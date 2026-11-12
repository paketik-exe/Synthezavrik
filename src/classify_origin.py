#!/usr/bin/env python3
"""Классификация эволюционного происхождения гена (алгоритм проекта «Синтезаврики»).

 1. Есть ли гомологи гена?            да → путь гомологии (2), нет → путь de novo (3)
 2. Где находятся гомологи?
    2.1 в том же геноме                          → Дупликация
    2.2 в далёких видах в синтенной позиции      → Древний ген
    2.3 в далёких видах вне синтении:
        2.3.1 мозаичная гомология                → Слияние / перестройка
        2.3.2 гомолог из другой ветви жизни      → HGT
              иначе                              → Неопределённо
 3. Есть ли синтенное некодирующее выравнивание в аутгруппе?
    3.1 да                                       → De novo
    3.2 связь с TE / ncRNA                       → Иное происхождение (TE / ncRNA)
        иначе                                    → Неопределённо

Вызывается из classify_gene.sh. Форматы таблиц (BLAST outfmt 6):
  белковые:     qseqid sseqid pident length qcovs evalue bitscore qlen slen qstart qend sstart send sseq
  нуклеотидные: те же без sseq
Источник/вид задаются именем файла:
  *_cand_prot.tsv, *_cand_nt.tsv        датасет кандидатов (вид = префикс ID)
  *_refprot_<вид>.tsv                   аннотированный протеом
  *_prot_refgen_<вид>.tsv, *_nt_refgen_<вид>.tsv   геном (tblastn / blastn)
  *_hgt_<группа>.tsv                    белки другой ветви жизни (бактерии, вирусы, ...)
  *_te.tsv, *_ncrna.tsv                 библиотеки мобильных элементов и ncRNA (blastn гена)
  synteny_<вид>.tsv                     синтенные окна у вида (из фланков гена)
"""
import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path

V = {
    "duplication": "Дупликация",
    "ancient": "Древний ген",
    "fusion": "Иные механизмы: слияние / перестройка",
    "hgt": "Иные механизмы: горизонтальный перенос (HGT)",
    "denovo": "De novo",
    "te_ncrna": "Иные механизмы: одомашнивание мобильного элемента / из ncRNA",
    "undetermined": "Неопределённо",
}


@dataclass
class Hit:
    source: str
    sseqid: str
    gene: str
    species: str
    pident: float
    length: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    sseq: str = ""
    syntenic: bool | None = None

    @property
    def qcov(self) -> float:
        return 100.0 * (self.qend - self.qstart + 1) / self.qlen if self.qlen else 0.0

    @property
    def scov(self) -> float:
        return 100.0 * (abs(self.send - self.sstart) + 1) / self.slen if self.slen else 0.0

    @property
    def stops(self) -> int:
        return self.sseq.count("*")

    @property
    def intact(self) -> bool:
        return self.source in ("refprot", "hgt") or self.stops == 0

    def short(self) -> str:
        syn = {True: "синтенно", False: "вне синтении", None: "синтения ?"}[self.syntenic]
        return (f"{self.sseqid} [{self.species}] {self.pident:.1f}% cov={self.qcov:.0f}% "
                f"e={self.evalue:.1e} bits={self.bitscore:.0f}"
                + (f" стопов={self.stops}" if self.sseq else "") + f" ({syn})")


# ------------------------------------------------------------------ чтение
def species_of(seq_id: str, rx) -> str:
    m = rx.search(seq_id)
    return m.group(1) if m else "?"


def read_hits(path: Path, rx, protein: bool) -> list[Hit]:
    stem = path.stem
    source, sp = "cand", None
    for tag, src in (("_refprot_", "refprot"), ("_prot_refgen_", "refgen"),
                     ("_nt_refgen_", "refgen"), ("_hgt_", "hgt")):
        if tag in stem:
            source, sp = src, stem.split(tag, 1)[1]
            break
    else:
        if stem.endswith("_te"):
            source, sp = "te", "TE"
        elif stem.endswith("_ncrna"):
            source, sp = "ncrna", "ncRNA"
    hits = []
    if not path.exists():
        return hits
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 13:
            continue
        gene = re.sub(r"_[1-6]$", "", f[1]) if (source == "cand" and protein) else f[1]
        q1, q2 = int(f[9]), int(f[10])
        hits.append(Hit(source, f[1], gene, sp or species_of(gene, rx), float(f[2]), int(f[3]),
                        float(f[5]), float(f[6]), int(f[7]), int(f[8]), min(q1, q2), max(q1, q2),
                        int(f[11]), int(f[12]), f[13] if protein and len(f) > 13 else ""))
    return hits


def read_windows(paths) -> dict[str, list[tuple[str, int, int, str]]]:
    win: dict[str, list] = {}
    for p in paths:
        p = Path(p)
        if not p.exists():
            continue
        for line in p.read_text().splitlines()[1:]:
            f = line.split("\t")
            if len(f) >= 5 and f[4] != "not_found":
                win.setdefault(f[0], []).append((f[1], int(f[2]), int(f[3]), f[4]))
    return win


# ------------------------------------------------------------------ дерево
@dataclass
class Node:
    name: str = ""
    children: list = None
    parent: "Node" = None

    def __post_init__(self):
        self.children = self.children or []

    def leaves(self):
        if not self.children:
            yield self
        for c in self.children:
            yield from c.leaves()


def parse_newick(text: str) -> Node:
    text = text.strip().rstrip(";")
    root = cur = Node()
    stack, tok = [], ""
    for ch in text:
        if ch == "(":
            child = Node(parent=cur); cur.children.append(child); stack.append(cur); cur = child
        elif ch == ",":
            cur.name = tok.split(":")[0].strip("' "); tok = ""
            child = Node(parent=stack[-1]); stack[-1].children.append(child); cur = child
        elif ch == ")":
            cur.name = tok.split(":")[0].strip("' "); tok = ""; cur = stack.pop()
        else:
            tok += ch
    cur.name = tok.split(":")[0].strip("' ")
    return root


def tree_sister(treefile: str, rec_id: str, rx):
    p = Path(treefile) if treefile else None
    if not p or not p.exists():
        return None
    root = parse_newick(p.read_text())
    q = next((l for l in root.leaves() if l.name.startswith(rec_id + "_frame")), None)
    if q is None or q.parent is None:
        return None
    par = q.parent
    sis = [l.name for c in par.children if c is not q for l in c.leaves()]
    nums = re.findall(r"[\d.]+", par.name) if par.parent is not None else []
    alrt, uf = (float(nums[0]), float(nums[1])) if len(nums) >= 2 else (None, float(nums[0]) if nums else None)
    sp = sorted({species_of(re.sub(r"_[1-6]_\d+-\d+$", "", s), rx) for s in sis})
    return sp, alrt, uf


# ------------------------------------------------------------------ мозаика
def mosaic(hits, part_cov, full_cov, max_overlap):
    best: dict[str, Hit] = {}
    for h in hits:
        if part_cov <= h.qcov < full_cov and (h.gene not in best or h.bitscore > best[h.gene].bitscore):
            best[h.gene] = h
    chosen = []
    for h in sorted(best.values(), key=lambda x: -x.bitscore):
        if all(100.0 * max(0, min(h.qend, c.qend) - max(h.qstart, c.qstart) + 1)
               / min(h.qend - h.qstart + 1, c.qend - c.qstart + 1) <= max_overlap for c in chosen):
            chosen.append(h)
    cov = set()
    for c in chosen:
        cov.update(range(c.qstart, c.qend + 1))
    return chosen, (100.0 * len(cov) / chosen[0].qlen if chosen else 0.0)


# ------------------------------------------------------------------ main
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rec-id", required=True)
    ap.add_argument("--focal", required=True)
    ap.add_argument("--ingroup", default="")
    ap.add_argument("--outgroup", default="")
    ap.add_argument("--species-regex", default=r"^([^_]+)_")
    ap.add_argument("--gene-key-regex", default=r"(FBgn\d+)")
    ap.add_argument("--orf-tsv", required=True)
    ap.add_argument("--prot-hits", nargs="*", default=[])
    ap.add_argument("--nt-hits", nargs="*", default=[])
    ap.add_argument("--hgt-hits", nargs="*", default=[])
    ap.add_argument("--te-hits", nargs="*", default=[])
    ap.add_argument("--ncrna-hits", nargs="*", default=[])
    ap.add_argument("--synteny", nargs="*", default=[])
    ap.add_argument("--orthologs", default="",
                    help="TSV ортологов NCBI (species, ..., protein): белок-ортолог = синтенный гомолог")
    ap.add_argument("--treefile", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--require-synteny", action="store_true",
                    help="не считать синтенными хиты, синтению которых проверить нельзя")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--min-pident", type=float, default=30)
    ap.add_argument("--min-cov", type=float, default=50)
    ap.add_argument("--min-orf-aa", type=int, default=30)
    ap.add_argument("--part-cov", type=float, default=20)
    ap.add_argument("--max-overlap", type=float, default=20)
    ap.add_argument("--nt-pident", type=float, default=70)
    ap.add_argument("--nt-cov", type=float, default=50)
    ap.add_argument("--min-ai", type=float, default=45)
    ap.add_argument("--te-pident", type=float, default=80)
    ap.add_argument("--te-evalue", type=float, default=1e-10)
    ap.add_argument("--te-cov", type=float, default=50)
    ap.add_argument("--self-pident", type=float, default=99,
                    help="хит в референсе фокального вида с такой идентичностью на ≥95%% длины = сам ген")
    ap.add_argument("--min-uf", type=float, default=95)
    ap.add_argument("--min-alrt", type=float, default=80)
    a = ap.parse_args()

    rx = re.compile(a.species_regex)
    krx = re.compile(a.gene_key_regex)
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    R: list[str] = [f"# Происхождение гена {a.rec_id}\n"]
    say = R.append

    rows = Path(a.orf_tsv).read_text().splitlines()
    orf = dict(zip(rows[0].split("\t"), rows[1].split("\t")))
    orf_aa = int(orf.get("orf_len_aa", 0))
    o1, o2 = sorted((int(orf.get("start", 0)), int(orf.get("end", 0))))

    key_q = (m.group(1) if (m := krx.search(a.rec_id)) else None)
    windows = read_windows(a.synteny)
    orth: dict[str, set] = {}
    if a.orthologs and Path(a.orthologs).exists():
        for line in Path(a.orthologs).read_text().splitlines()[1:]:
            f = line.split("\t")
            if len(f) >= 5:
                orth.setdefault(f[0], set()).add(f[4])
    acc_rx = re.compile(r"[A-Z]{2}_\d+\.\d+")
    acc_of = lambda sid: (m.group(0) if (m := acc_rx.search(sid)) else sid)

    dup_records: list[str] = []

    def is_self(h):
        if h.source == "cand":
            if h.gene == a.rec_id:
                return True
            # другая запись датасета того же вида, идентичная почти на всю длину, —
            # это тот же участок генома (дубль/перекрытие записей), а не паралог
            if species_of(h.gene, rx) == a.focal and h.pident >= a.self_pident and h.qcov >= 95:
                dup_records.append(h.gene)
                return True
            return False
        # в протеоме/геноме фокального вида есть и сам ген: узнаём его по ключу
        # (FBgn в ID субъекта) или по почти полной идентичности на всю длину
        if h.species == a.focal and h.source in ("refprot", "refgen"):
            if h.source == "refprot" and a.focal in orth:      # белки самого гена из NCBI
                return acc_of(h.sseqid) in orth[a.focal]
            k = krx.search(h.sseqid)
            if k and key_q:
                return k.group(1) == key_q
            return h.pident >= a.self_pident and h.qcov >= 95
        return False

    def synteny_of(h: Hit):
        if h.source == "refgen" and h.species in windows:
            s1, s2 = sorted((h.sstart, h.send))
            return any(c == h.sseqid and s1 <= w2 and s2 >= w1 for c, w1, w2, _ in windows[h.species])
        if h.source == "refprot" and h.species in orth:       # ортолог NCBI = синтения
            return acc_of(h.sseqid) in orth[h.species]
        if h.source == "cand" and key_q:
            k = krx.search(h.gene)
            return bool(k) and k.group(1) == key_q
        return None

    prot = [h for p in a.prot_hits for h in read_hits(Path(p), rx, True) if not is_self(h)]
    nt = [h for p in a.nt_hits for h in read_hits(Path(p), rx, False) if not is_self(h)]
    hgt = [h for p in a.hgt_hits for h in read_hits(Path(p), rx, True)]
    te = [h for p in a.te_hits for h in read_hits(Path(p), rx, False)]
    nc = [h for p in a.ncrna_hits for h in read_hits(Path(p), rx, False)]
    for h in prot + nt:
        h.syntenic = synteny_of(h)
    if dup_records:
        say(f"Записи датасета, совпадающие с геном на ≥ {a.self_pident:g}% по ≥ 95% длины "
            f"(считаются тем же участком, не паралогом): {', '.join(sorted(set(dup_records)))}\n")
    syn_ok = (lambda h: h.syntenic is True) if a.require_synteny else (lambda h: h.syntenic is not False)

    all_sp = {h.species for h in prot + nt} - {a.focal}
    ingroup = {s for s in a.ingroup.split(",") if s}
    outgroup = {s for s in a.outgroup.split(",") if s} or (all_sp - ingroup)
    best = lambda hs: max(hs, key=lambda h: h.bitscore) if hs else None

    # значимые белковые гомологи с целой рамкой (без ingroup-ортологов самого гена)
    sig = [h for h in prot if h.evalue <= a.evalue and h.pident >= a.min_pident and h.intact]
    sig_hgt = [h for h in hgt if h.evalue <= a.evalue and h.pident >= a.min_pident]
    homologs = [h for h in sig if h.qcov >= a.part_cov and h.species not in ingroup] + \
               [h for h in sig_hgt if h.qcov >= a.part_cov]
    full = [h for h in sig if h.qcov >= a.min_cov]
    ing_intact = sorted({h.species for h in full if h.species in ingroup})

    m: dict = {"gene": a.rec_id, "focal": a.focal, "orf_len_aa": orf_aa,
               "n_prot_hits": len(prot), "n_homologs": len(homologs),
               "ingroup_intact": ",".join(ing_intact) or "-"}
    say(f"Фокальный вид **{a.focal}**; ingroup: {', '.join(sorted(ingroup)) or '—'}; "
        f"outgroup (далёкие): {', '.join(sorted(outgroup)) or '—'}; "
        f"синтения: окна для {', '.join(sorted(windows)) or '—'}, ортологи NCBI у "
        f"{', '.join(sorted(set(orth) - {a.focal})) or '—'}, ключ гена {key_q or '—'}\n")
    say(f"Пороги: e ≤ {a.evalue:g}, pident ≥ {a.min_pident:g}%, полноразмерный ≥ {a.min_cov:g}%, "
        f"часть ≥ {a.part_cov:g}%; ДНК ≥ {a.nt_pident:g}% / ≥ {a.nt_cov:g}%; AI ≥ {a.min_ai:g}; "
        f"TE/ncRNA: e ≤ {a.te_evalue:g}, ≥ {a.te_pident:g}%, перекрытие ≥ {a.te_cov:g}%\n")

    verdict = conf = reason = None

    # ---------------- 0. ORF
    say("## 0. Белок-кодирующая рамка")
    say(f"- самый длинный ORF: **{orf_aa} а.к.** ({orf.get('orf_id', '—')}, {orf.get('orf_frac', '?')}% длины)")
    if orf_aa < a.min_orf_aa:
        verdict, conf, reason = "undetermined", "—", f"нет ORF ≥ {a.min_orf_aa} а.к. — не белок-кодирующий"
        say(f"- ✗ {reason}\n")

    # ---------------- 1. Есть ли гомологи
    if verdict is None:
        say("\n## 1. Есть ли гомологи гена?")
        say(f"- значимых белковых гомологов с целой рамкой (вне ingroup, покрытие ≥ {a.part_cov:g}%): "
            f"**{len(homologs)}** (в т.ч. другой ветви жизни: {len([h for h in homologs if h.source == 'hgt'])})")
        if ing_intact:
            say(f"- ортологи с целой рамкой у ingroup (не учитываются — возраст гена): {', '.join(ing_intact)}")
        path = 2 if homologs else 3
        say(f"- → {'1.1 да: не de novo, путь гомологии' if path == 2 else '1.2 нет: кандидат de novo'}")

        if path == 2:
            say("\n## 2. Где находятся гомологи?")
            par = [h for h in full if h.species == a.focal]
            anc = [h for h in full if h.species in outgroup and syn_ok(h)]
            far_nonsyn = [h for h in homologs if h.source != "hgt" and h.species in outgroup and h.syntenic is False]
            bp, ba = best(par), best(anc)
            m.update(paralogs=len({h.gene for h in par}), par_best=bp.short() if bp else "-",
                     anc_species=",".join(sorted({h.species for h in anc})) or "-",
                     anc_best=ba.short() if ba else "-")
            say(f"- 2.1 паралоги в геноме {a.focal}: {m['paralogs']}" + (f"; лучший: {bp.short()}" if bp else ""))
            say(f"- 2.2 полноразмерные гомологи у далёких видов в синтенной позиции: {m['anc_species']}"
                + (f"; лучший: {ba.short()}" if ba else ""))
            ts = tree_sister(a.treefile, a.rec_id, rx)
            if ts:
                say(f"- дерево: сестринская группа — {', '.join(ts[0])}; SH-aLRT={ts[1]}, UFBoot={ts[2]}")
            supp = bool(ts) and ((ts[2] or 0) >= a.min_uf or (ts[1] or 0) >= a.min_alrt)

            if par and anc:
                # Уточнение п. 2.1/2.2: древний ген может иметь древний паралог.
                if supp and set(ts[0]) & outgroup:
                    verdict, conf = "ancient", "высокая"
                    reason = "в дереве ген надёжно группируется с ортологом далёкого вида; паралог — более древняя дупликация"
                elif (supp and ts[0] == [a.focal]) or bp.bitscore > ba.bitscore:
                    verdict, conf = "duplication", "высокая" if supp else "средняя"
                    reason = (f"паралог ближе к гену (bits {bp.bitscore:.0f}), чем синтенный гомолог далёкого вида "
                              f"(bits {ba.bitscore:.0f}) → дупликация после расхождения с аутгруппой")
                else:
                    verdict, conf = "ancient", "средняя"
                    reason = (f"синтенный гомолог далёкого вида ближе (bits {ba.bitscore:.0f}), чем паралог "
                              f"(bits {bp.bitscore:.0f}) → ген старше дупликации")
            elif par:
                verdict, conf = "duplication", "высокая" if supp else "средняя"
                reason = f"2.1 полноразмерный паралог в том же геноме: {bp.short()}"
            elif anc:
                n_syn_true = len({h.species for h in anc if h.syntenic is True})
                verdict = "ancient"
                conf = "высокая" if n_syn_true >= 1 else "средняя (синтения не проверена)"
                reason = f"2.2 гомологи с целой рамкой у далёких видов в синтенной позиции: {m['anc_species']}"
            if verdict:
                say(f"- ✓ {V[verdict]}")
            else:
                say(f"- 2.3 гомологи только вне синтении / частичные (вне синтении: {len(far_nonsyn)})")
                parts, union = mosaic([h for h in homologs if h.source != "hgt"],
                                      a.part_cov, a.min_cov, a.max_overlap)
                m.update(mosaic_parts=len(parts), mosaic_cov=round(union, 1))
                say(f"- 2.3.1 мозаичная гомология: частей {len(parts)}, суммарное покрытие {union:.0f}%")
                for p in parts:
                    say(f"    - {p.qstart}–{p.qend}: {p.short()}")
                if len(parts) >= 2 and union >= a.min_cov:
                    verdict, conf = "fusion", "средняя"
                    reason = f"{len(parts)} разных гена гомологичны разным участкам гена (суммарно {union:.0f}%)"
                else:
                    own = min([h.evalue for h in prot if h.evalue <= 10] or [10.0])
                    oth = best(sig_hgt)
                    ai = (math.log(own + 1e-200) - math.log(min(h.evalue for h in sig_hgt) + 1e-200)
                          if sig_hgt else None)
                    m.update(hgt_best=oth.short() if oth else "-", alien_index=round(ai, 1) if ai is not None else "-")
                    say(f"- 2.3.2 гомолог из другой ветви жизни: {oth.short() if oth else 'нет'}; "
                        f"Alien Index = {('%.1f' % ai) if ai is not None else '—'}")
                    if oth and ai is not None and ai > 0:
                        verdict = "hgt"
                        conf = "высокая" if ai >= a.min_ai else "низкая"
                        reason = (f"лучший гомолог — в другой ветви жизни ({oth.species}), он сильнее любого "
                                  f"гомолога своей ветви: AI = {ai:.1f}")
                    else:
                        verdict, conf = "undetermined", "—"
                        reason = ("гомологи есть, но не в своём геноме, не в синтенной позиции у далёких видов, "
                                  "не мозаичны и не указывают на HGT")
                say(f"- → {V[verdict]}")

        # ---------------- 3. Путь de novo
        if path == 3:
            say("\n## 3. Есть ли синтенное некодирующее выравнивание в аутгруппе?")
            nt_og = [h for h in nt if h.species in outgroup and h.evalue <= a.evalue
                     and h.pident >= a.nt_pident and h.qcov >= a.nt_cov and syn_ok(h)]
            noncod = []
            for sp in sorted({h.species for h in nt_og}):
                bn = best([h for h in nt_og if h.species == sp])
                bpr = best([h for h in prot if h.species == sp and h.syntenic is not False])
                coding = bpr is not None and bpr.intact and bpr.qcov >= a.min_cov and bpr.evalue <= a.evalue
                state = ("белкового сходства нет" if bpr is None else
                         f"стопов {bpr.stops}, покрытие {bpr.qcov:.0f}%")
                say(f"- {sp}: ДНК {bn.pident:.1f}% / {bn.qcov:.0f}% ({'синтенно' if bn.syntenic else 'синтения ?'}); "
                    f"рамка: {state} → {'кодирующий' if coding else 'некодирующий'}")
                if not coding:
                    noncod.append((sp, bn))
            m.update(noncoding_outgroup=",".join(s for s, _ in noncod) or "-")
            if noncod:
                verdict = "denovo"
                checked = any(bn.syntenic is True for _, bn in noncod)
                conf = ("высокая" if checked and len(noncod) >= 2 else
                        "средняя" if checked else "низкая (синтения не проверена)")
                reason = (f"3.1 у аутгруппы ({', '.join(s for s, _ in noncod)}) есть синтенный гомологичный "
                          f"участок ДНК без целой рамки"
                          + (f"; рамка цела у ingroup: {', '.join(ing_intact)}" if ing_intact else ""))
                say(f"- ✓ {reason}")
            else:
                say("- ✗ синтенного некодирующего выравнивания нет")
                say("\n## 3.2 Есть ли связь с TE / ncRNA?")

                def link(hs, label):
                    good = []
                    for h in hs:
                        ov = max(0, min(h.qend, o2) - max(h.qstart, o1) + 1) if o2 else 0
                        orf_ov = 100.0 * ov / (o2 - o1 + 1) if o2 > o1 else 0.0
                        if h.evalue <= a.te_evalue and h.pident >= a.te_pident and \
                                (orf_ov >= a.te_cov or h.scov >= a.te_cov):
                            good.append((h, orf_ov))
                    for h, ov in sorted(good, key=lambda x: -x[0].bitscore)[:3]:
                        say(f"    - {label}: {h.sseqid} {h.pident:.1f}%, e={h.evalue:.1e}, "
                            f"перекрытие с ORF {ov:.0f}%, покрытие элемента {h.scov:.0f}%")
                    return good

                te_ok, nc_ok = link(te, "TE"), link(nc, "ncRNA")
                m.update(te_link=len(te_ok), ncrna_link=len(nc_ok))
                say(f"- TE: {len(te_ok)}, ncRNA: {len(nc_ok)}"
                    + ("" if (te or nc or a.te_hits or a.ncrna_hits) else " (библиотеки TE/ncRNA не заданы)"))
                if te_ok or nc_ok:
                    verdict, conf = "te_ncrna", "средняя"
                    reason = ("3.2.1 ORF перекрывается с " + " и ".join(
                        x for x, ok in (("мобильным элементом", te_ok), ("ncRNA", nc_ok)) if ok))
                else:
                    verdict, conf = "undetermined", "—"
                    reason = ("3.2.2 нет ни белковых гомологов, ни синтенного некодирующего участка у аутгруппы, "
                              "ни связи с TE/ncRNA" + ("" if outgroup else "; не задано ни одного вида аутгруппы"))
                say(f"- → {V[verdict]}")

    m.update(verdict=V[verdict], confidence=conf, reason=reason)
    say(f"\n## Итог: **{V[verdict]}** (уверенность: {conf})\n\n{reason}\n")
    keys = list(m)
    (out / "origin.tsv").write_text("\t".join(keys) + "\n" + "\t".join(str(m[k]) for k in keys) + "\n")
    (out / "origin_report.md").write_text("\n".join(R))
    print(f"{a.rec_id}\t{V[verdict]}\t{conf}\t{reason}")


if __name__ == "__main__":
    main()
