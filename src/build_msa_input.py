#!/usr/bin/env python3
"""Собирает вход для MAFFT только из выровненных BLAST участков.

Использование:
  build_msa_input.py <top.tsv> <6frames.faa> <лучшая_рамка> <out.faa> <hit_ids.txt>

top.tsv — вывод blastp в формате
  qseqid sseqid pident length qcovs evalue bitscore stitle qstart qend sstart send sseq
(уже отфильтрованный, отсортированный по bitscore и только для лучшей рамки).

Берём:
  * у хитов — только выровненный участок (sseq без гэпов), а не весь белок:
    иначе длинная геномная последовательность даёт «белок» в миллионы а.к.
    и MAFFT падает по памяти;
  * у запроса — участок лучшей рамки от min(qstart) до max(qend) по хитам.
Стоп-кодоны '*' заменяются на 'X' (MAFFT/IQ-TREE их не любят).
"""
import sys

from Bio import SeqIO


def main() -> None:
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    top_tsv, frames_faa, best_frame, out_faa, ids_out = sys.argv[1:]

    hits, qs, qe = [], None, None
    with open(top_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                sys.exit(f"ОШИБКА: в {top_tsv} нет колонок qstart/qend/sstart/send/sseq")
            sseqid = f[1]
            q1, q2, s1, s2 = (int(x) for x in f[8:12])
            seq = f[12].replace("-", "").replace("*", "X")
            hits.append((f"{sseqid}_{s1}-{s2}", sseqid, seq))
            qs = q1 if qs is None else min(qs, q1)
            qe = q2 if qe is None else max(qe, q2)

    if not hits:
        sys.exit("ОШИБКА: нет хитов для лучшей рамки")

    query = next((r for r in SeqIO.parse(frames_faa, "fasta") if r.id == best_frame), None)
    if query is None:
        sys.exit(f"ОШИБКА: рамка {best_frame} не найдена в {frames_faa}")
    qseq = str(query.seq[qs - 1:qe]).replace("*", "X")

    with open(out_faa, "w") as fh:
        fh.write(f">{best_frame}_{qs}-{qe}\n{qseq}\n")
        for rid, _, seq in hits:
            fh.write(f">{rid}\n{seq}\n")
    with open(ids_out, "w") as fh:
        fh.write("\n".join(sid for _, sid, _ in hits) + "\n")

    lens = [len(s) for *_, s in hits]
    print(f"Запрос: {best_frame} [{qs}-{qe}], {len(qseq)} а.к.; "
          f"хитов: {len(hits)}, длины {min(lens)}–{max(lens)} а.к.", file=sys.stderr)


if __name__ == "__main__":
    main()
