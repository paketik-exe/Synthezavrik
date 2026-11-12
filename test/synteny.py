#!/usr/bin/env python3
"""Синтения через фланки гена.

  synteny.py locus <focal.gff> <ключ_гена>
      Ищет в GFF ген (feature=gene), в атрибутах которого есть ключ (например, FBgn0030912).
      Печатает: chrom<TAB>start<TAB>end<TAB>strand

  synteny.py window <flanks_vs_genome.tsv> <вид> <длина_гена> <фланк> <out.tsv>
      flanks_vs_genome.tsv — blastn фланков (ID left/right) против генома вида, outfmt:
        6 qseqid sseqid pident length evalue bitscore sstart send
      Находит пару хитов левого и правого фланка на одном контиге на разумном расстоянии.
      Окно между ними — синтенная позиция гена у этого вида.
      out.tsv: species contig start end status   (status: both | left | right | not_found)
"""
import re
import sys


def locus(gff: str, key: str) -> None:
    rx = re.compile(rf"(^|[=:,;\s]){re.escape(key)}([;,\s]|$)")
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[2] == "gene" and rx.search(f[8]):
                print(f"{f[0]}\t{f[3]}\t{f[4]}\t{f[6]}")
                return
    sys.exit(f"ОШИБКА: ген {key} не найден в {gff}")


def window(tsv: str, species: str, gene_len: int, flank: int, out: str) -> None:
    hits = {"left": [], "right": []}
    with open(tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 8 and f[0] in hits:
                s1, s2 = sorted((int(f[6]), int(f[7])))
                hits[f[0]].append((f[1], s1, s2, float(f[5])))
    max_span = 3 * (2 * flank + gene_len)
    best = None
    for c1, a1, b1, x1 in hits["left"]:
        for c2, a2, b2, x2 in hits["right"]:
            if c1 != c2:
                continue
            lo, hi = min(a1, a2), max(b1, b2)
            if hi - lo <= max_span and (best is None or x1 + x2 > best[3]):
                best = (c1, lo, hi, x1 + x2)
    if best:
        row = (best[0], best[1], best[2], "both")
    else:
        one = [(side, h) for side in ("left", "right") for h in hits[side]]
        if one:
            side, (c, a, b, _) = max(one, key=lambda x: x[1][3])
            pad = gene_len + flank
            row = (c, max(1, a - pad), b + pad, side)
        else:
            row = ("-", 0, 0, "not_found")
    with open(out, "w") as fh:
        fh.write("species\tcontig\tstart\tend\tstatus\n")
        fh.write(f"{species}\t{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\n")
    print(f"{species}: {row[3]} {row[0]}:{row[1]}-{row[2]}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "locus":
        locus(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 7 and sys.argv[1] == "window":
        window(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6])
    else:
        sys.exit(__doc__)
