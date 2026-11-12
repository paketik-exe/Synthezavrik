#!/usr/bin/env python3
"""Онлайн-доступ к NCBI для classify_gene.sh --online (только стандартная библиотека Python).

Подкоманды:
  geneid  <ключ> <организм>                 → GeneID (E-utilities esearch, база gene)
  locus   <GeneID>                          → acc<TAB>start<TAB>end  (координаты гена, 1-based)
  flanks  <acc> <start> <end> <N> <out.fna> → левый/правый фланки по N нт (efetch nuccore)
  orthologs <GeneID> <out.tsv>              → ортологи NCBI (сходство + синтения) и их белки
  split   <blast.tsv> <prefix> <kind> [--focal ВИД]
          раскладывает результат удалённого BLAST по видам в файлы формата classify_origin.py.
          Последняя колонка входа — stitle; вид берётся из [Genus species] или начала заголовка.
          kind: refprot | prot_refgen | nt_refgen | hgt | ncrna | flanks
  taxa    <вид,вид,...>                     → Entrez-фильтр "(Drosophila yakuba[Organism] OR ...)"
  name    <вид>                             → латинское название (dmel → Drosophila melanogaster)

Ключ NCBI (необязательно, ускоряет): export NCBI_API_KEY=...
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
TAXON = {
    "dmel": "Drosophila melanogaster", "dsim": "Drosophila simulans", "dsec": "Drosophila sechellia",
    "dmau": "Drosophila mauritiana", "dyak": "Drosophila yakuba", "dere": "Drosophila erecta",
    "dtei": "Drosophila teissieri", "dsan": "Drosophila santomea", "dana": "Drosophila ananassae",
    "dpse": "Drosophila pseudoobscura", "dper": "Drosophila persimilis", "dwil": "Drosophila willistoni",
    "dmoj": "Drosophila mojavensis", "dvir": "Drosophila virilis", "dgri": "Drosophila grimshawi",
    "dbus": "Drosophila busckii", "dsuz": "Drosophila suzukii", "dbia": "Drosophila biarmipes",
}
ACC_RE = re.compile(r"[A-Z]{2}_\d+\.\d+")


def get(url: str, params: dict | None = None, tries: int = 4) -> str:
    params = dict(params or {})
    if os.environ.get("NCBI_API_KEY") and "eutils" in url:
        params["api_key"] = os.environ["NCBI_API_KEY"]
    full = url + ("?" + urllib.parse.urlencode(params) if params else "")
    for i in range(tries):
        try:
            req = urllib.request.Request(full, headers={"Accept": "application/json",
                                                        "User-Agent": "synthezavrik/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                time.sleep(0.4)          # E-utilities: не более ~3 запросов в секунду
                return r.read().decode("utf-8", "replace")
        except Exception as e:           # noqa: BLE001 — сеть: повторяем
            if i == tries - 1:
                sys.exit(f"ОШИБКА сети: {e} ({full[:150]})")
            time.sleep(3 * (i + 1))
    return ""


def abbr(scientific: str) -> str:
    """'Drosophila yakuba' → 'dyak' (соглашение FlyBase)."""
    for k, v in TAXON.items():
        if v.lower() == scientific.lower():
            return k
    w = scientific.split()
    return (w[0][0] + w[1][:3]).lower() if len(w) >= 2 else scientific.lower()


def species_from_title(title: str) -> str:
    m = re.findall(r"\[([A-Z][a-z]+ [a-z]+)[^\]]*\]", title)
    if m:
        return m[-1]
    m = re.match(r"\s*([A-Z][a-z]+ [a-z]+)", title)
    return m.group(1) if m else "unknown species"


# ------------------------------------------------------------------ подкоманды
def cmd_geneid(key: str, organism: str) -> None:
    js = json.loads(get(f"{EUTILS}/esearch.fcgi", {"db": "gene", "retmode": "json",
                                                   "term": f'{key}[All Fields] AND "{organism}"[Organism]'}))
    ids = js.get("esearchresult", {}).get("idlist", [])
    if not ids:
        sys.exit(f"ОШИБКА: {key} не найден в NCBI Gene")
    print(ids[0])


def cmd_locus(gene_id: str) -> None:
    js = json.loads(get(f"{EUTILS}/esummary.fcgi", {"db": "gene", "id": gene_id, "retmode": "json"}))
    info = js["result"][gene_id].get("genomicinfo") or []
    if not info:
        sys.exit(f"ОШИБКА: у GeneID {gene_id} нет координат")
    g = info[0]
    s, e = sorted((int(g["chrstart"]) + 1, int(g["chrstop"]) + 1))
    print(f"{g['chraccver']}\t{s}\t{e}")


def cmd_flanks(acc: str, start: int, end: int, n: int, out: str) -> None:
    def seq(a, b):
        txt = get(f"{EUTILS}/efetch.fcgi", {"db": "nuccore", "id": acc, "rettype": "fasta",
                                             "retmode": "text", "seq_start": a, "seq_stop": b})
        return "".join(l for l in txt.splitlines() if not l.startswith(">"))
    left, right = seq(max(1, start - n), max(1, start - 1)), seq(end + 1, end + n)
    with open(out, "w") as fh:
        fh.write(f">left\n{left}\n>right\n{right}\n")
    print(f"фланки: {len(left)} + {len(right)} нт", file=sys.stderr)


def cmd_orthologs(gene_id: str, out: str) -> None:
    genes, token = [], None
    while True:
        p = {"page_size": 1000}
        if token:
            p["page_token"] = token
        js = json.loads(get(f"{DATASETS}/gene/id/{gene_id}/orthologs", p))
        for r in js.get("reports", []):
            g = r.get("gene", r)
            genes.append((str(g.get("gene_id")), g.get("taxname", "?"), g.get("symbol", "")))
        token = js.get("next_page_token")
        if not token:
            break
    rows = []
    for i in range(0, len(genes), 100):
        chunk = genes[i:i + 100]
        info = {gid: (tax, sym) for gid, tax, sym in chunk}
        js = json.loads(get(f"{DATASETS}/gene/id/{','.join(g for g, _, _ in chunk)}/product_report"))
        for r in js.get("reports", []):
            prod = r.get("product", r)
            gid = str(prod.get("gene_id"))
            tax, sym = info.get(gid, (prod.get("taxname", "?"), prod.get("symbol", "")))
            accs = sorted({a for a in ACC_RE.findall(json.dumps(prod)) if a[1] == "P"})
            for a in accs or ["-"]:
                rows.append((abbr(tax), tax, gid, sym, a))
    with open(out, "w") as fh:
        fh.write("species\ttaxname\tgene_id\tsymbol\tprotein\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"ортологов NCBI: {len(genes)} генов у {len({r[0] for r in rows})} видов", file=sys.stderr)


def cmd_split(path: str, prefix: str, kind: str, focal: str = "") -> None:
    groups: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            sp = abbr(species_from_title(f[-1]))
            m = ACC_RE.search(f[1])
            f[1] = m.group(0) if m else f[1]
            groups.setdefault(sp, []).append("\t".join(f[:-1]))
    for sp, rows in groups.items():
        name = {"refprot": f"{prefix}_refprot_{sp}.tsv", "prot_refgen": f"{prefix}_prot_refgen_{sp}.tsv",
                "nt_refgen": f"{prefix}_nt_refgen_{sp}.tsv", "flanks": f"{prefix}_flanks_{sp}.tsv",
                "hgt": f"{prefix}_hgt_other_life.tsv", "ncrna": f"{prefix}_ncrna.tsv"}[kind]
        with open(name, "a") as fh:
            fh.write("\n".join(rows) + "\n")
    print(f"{kind}: хиты у {len(groups)} видов ({', '.join(sorted(groups)) or '—'})", file=sys.stderr)


def cmd_taxa(spec: str) -> None:
    names = []
    for s in filter(None, spec.split(",")):
        if s not in TAXON:
            sys.exit(f"ОШИБКА: неизвестный вид '{s}'. Известны: {', '.join(TAXON)}")
        names.append(f'"{TAXON[s]}"[Organism]')
    print("(" + " OR ".join(names) + ")")


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a:
        sys.exit(__doc__)
    c = a[0]
    if c == "geneid" and len(a) == 3:
        cmd_geneid(a[1], a[2])
    elif c == "locus" and len(a) == 2:
        cmd_locus(a[1])
    elif c == "flanks" and len(a) == 6:
        cmd_flanks(a[1], int(a[2]), int(a[3]), int(a[4]), a[5])
    elif c == "orthologs" and len(a) == 3:
        cmd_orthologs(a[1], a[2])
    elif c == "split" and len(a) in (4, 6):
        cmd_split(a[1], a[2], a[3], a[5] if len(a) == 6 else "")
    elif c == "name" and len(a) == 2:
        print(TAXON.get(a[1]) or sys.exit(f"ОШИБКА: неизвестный вид '{a[1]}'"))
    elif c == "taxa" and len(a) == 2:
        cmd_taxa(a[1])
    else:
        sys.exit(__doc__)
