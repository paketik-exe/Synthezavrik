"""Веб-интерфейс пайплайна: загрузил файл -> выбрал ген -> нажал кнопку.

Запуск:  ./run_app.sh   (затем открыть http://localhost:8501)
"""
import importlib.util
import io
import os
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts" if (ROOT / "scripts" / "txt_to_fasta.py").exists() else ROOT
UPLOADS = ROOT / "data" / "uploads"
RESULTS = ROOT / "results"
BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "qcovs", "evalue", "bitscore", "stitle"]
REQUIRED = ["blastn", "blastp", "makeblastdb", "blastdbcmd", "mafft", "transeq", "getorf"]


def load_parser():
    spec = importlib.util.spec_from_file_location("txt_to_fasta", SCRIPTS / "txt_to_fasta.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.parse_records


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", s)


def zip_dir(path: Path) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in path.rglob("*"):
            if f.is_file():
                zf.write(f, f.relative_to(path.parent))
    return buf.getvalue()


def read_blast(path: Path) -> pd.DataFrame | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    return pd.read_csv(path, sep="\t", names=BLAST_COLS)


def show_df(df: pd.DataFrame | None) -> None:
    if df is None:
        st.write("Нет данных")
    else:
        st.dataframe(df, use_container_width=True)


def show_gene(dataset: str, gene: str) -> None:
    gs = safe_name(gene)
    d = RESULTS / dataset / gs
    if not d.exists():
        st.warning("Результатов нет.")
        return
    phylo = d / "phylogeny"

    best = d / f"{gs}_best_frame.txt"
    c1, c2, c3 = st.columns(3)
    c1.metric("Лучшая рамка", best.read_text().strip() if best.exists() else "—")
    ids = d / f"{gs}_blastp_hit_ids.txt"
    c2.metric("Хитов в дереве", len(ids.read_text().split()) if ids.exists() else 0)
    iq = phylo / f"{gs}_blastp.iqtree"
    model = "—"
    if iq.exists():
        m = re.search(r"Best-fit model.*?:\s*(\S+)", iq.read_text())
        model = m.group(1) if m else "—"
    c3.metric("Модель (IQ-TREE)", model)

    png = phylo / f"{gs}_blastp.png"
    if png.exists():
        st.image(str(png), caption="Филогенетическое дерево (поддержки: SH-aLRT / UFBoot)")
    else:
        st.info("Дерево не построено (нужно ≥ 4 последовательностей в выравнивании).")

    with st.expander("Топ хитов blastp (после фильтрации)"):
        df = next((read_blast(p) for p in sorted(d.glob(f"{gs}_blastp.top*.tsv"))), None)
        show_df(df)
    with st.expander("Все хиты blastn"):
        show_df(read_blast(d / f"{gs}_blastn.tsv"))
    ascii_tree = phylo / f"{gs}_blastp.ascii.txt"
    if ascii_tree.exists():
        with st.expander("Дерево в текстовом виде"):
            st.code(ascii_tree.read_text())
    orfs = d / f"{gs}_orfs.fasta"
    if orfs.exists():
        with st.expander(f"ORF ({orfs.read_text().count('>')} шт.)"):
            st.code(orfs.read_text()[:20000])

    st.download_button(f"Скачать все файлы по {gene} (zip)", zip_dir(d),
                       file_name=f"{dataset}_{gs}.zip", key=f"dl_{gs}")


# ============================== UI ==============================
st.set_page_config(page_title="Анализ гена", layout="wide")
st.title("Анализ гена: BLAST → выравнивание → филогения")

missing = [c for c in REQUIRED if shutil.which(c) is None]
if not (shutil.which("iqtree2") or shutil.which("iqtree")):
    missing.append("iqtree2")
if missing:
    st.error("Не найдены программы: " + ", ".join(missing) + ". Выполните `./setup.sh`.")
    st.stop()

with st.sidebar:
    st.header("Параметры")
    evalue = st.text_input("e-value", "1e-5")
    threads = st.number_input("Потоки", 1, os.cpu_count() or 1, min(4, os.cpu_count() or 1))
    min_pident = st.number_input("Мин. % идентичности", 0, 100, 30)
    min_qcov = st.number_input("Мин. покрытие запроса, %", 0, 100, 50)
    min_len = st.number_input("Мин. длина выравнивания", 1, 100000, 30)
    top = st.number_input("Сколько лучших хитов брать", 1, 1000, 50)
    max_targets = st.number_input("max_target_seqs", 1, 10000, 100)
    min_orf = st.number_input("Мин. длина ORF", 30, 100000, 90)

up = st.file_uploader("Файл с последовательностями (.txt с блоками «raw sequence for …» или FASTA)",
                      type=["txt", "fasta", "fa", "fna", "fas"])
if up is None:
    st.stop()

records = load_parser()(up.getvalue().decode("utf-8", errors="replace"))
if not records:
    st.error("В файле не найдено ни одной последовательности.")
    st.stop()

UPLOADS.mkdir(parents=True, exist_ok=True)
in_path = UPLOADS / safe_name(up.name)
in_path.write_bytes(up.getvalue())
dataset = safe_name(Path(up.name).stem)
ids = [rid for rid, _ in records]
st.success(f"Найдено последовательностей: {len(ids)} · датасет: `{dataset}`")

mode = st.radio("Что анализировать", ["Один ген", "Несколько генов", "Все гены файла"], horizontal=True)
if mode == "Один ген":
    genes = [st.selectbox("Ген", ids)]
elif mode == "Несколько генов":
    genes = st.multiselect("Гены", ids)
else:
    genes = ids
    st.caption(f"Будет обработано {len(ids)} генов — это может занять много времени.")

if st.button("▶ Запустить анализ", type="primary", disabled=not genes):
    cmd = ["bash", str(ROOT / "run_all.sh"), str(in_path), "-n", dataset,
           "-e", evalue, "-t", str(threads), "--min-pident", str(min_pident),
           "--min-qcov", str(min_qcov), "--min-len", str(min_len), "--top", str(top),
           "--max-targets", str(max_targets), "--min-orf", str(min_orf)]
    if mode == "Все гены файла":
        cmd.append("--all")
    elif len(genes) == 1:
        cmd += ["-g", genes[0]]
    else:
        gfile = UPLOADS / f"{dataset}_genes.txt"
        gfile.write_text("\n".join(genes) + "\n")
        cmd += ["--genes", str(gfile)]

    progress = st.progress(0.0, text="Запуск…")
    log_box = st.empty()
    lines: list[str] = []
    n = max(len(genes), 1)
    done_genes = 0
    proc = subprocess.Popen(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    for line in proc.stdout:
        lines.append(line.rstrip("\n"))
        if m := re.search(r"=== \[(\d+)/(\d+)\]", line):
            done_genes = int(m.group(1)) - 1
        if m := re.search(r"\[(\d+)/10\]", line):
            frac = (done_genes + int(m.group(1)) / 10) / n
            progress.progress(min(frac, 1.0), text=line.strip()[-120:])
        log_box.code("\n".join(lines[-40:]), language="text")
    rc = proc.wait()
    progress.progress(1.0, text="Готово" if rc == 0 else "Завершено с ошибками")
    st.session_state["last"] = {"dataset": dataset, "genes": genes, "rc": rc, "log": lines}

last = st.session_state.get("last")
if last and last["dataset"] == dataset:
    if last["rc"] == 0:
        st.success("Анализ завершён.")
    else:
        st.warning("Часть генов завершилась с ошибкой — см. лог.")
    with st.expander("Полный лог"):
        st.code("\n".join(last["log"]), language="text")

    summary = RESULTS / dataset / "summary.tsv"
    if summary.exists() and len(last["genes"]) > 1:
        st.subheader("Сводка")
        st.dataframe(pd.read_csv(summary, sep="\t"), use_container_width=True)

    shown = last["genes"] if len(last["genes"]) <= 20 else last["genes"][:20]
    tabs = st.tabs(shown)
    for tab, g in zip(tabs, shown):
        with tab:
            show_gene(dataset, g)

    archive = RESULTS / f"{dataset}.zip"
    if archive.exists():
        st.download_button("Скачать все результаты датасета (zip)", archive.read_bytes(),
                           file_name=archive.name)
