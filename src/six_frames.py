#!/usr/bin/env python3
"""Транслирует нуклеотидные последовательности в 6 рамках считывания.

Использование: six_frames.py <input.fasta> <output.faa>
ID белков: <ID>_frame<1-3><+|->  (без символа '|', чтобы BLAST не искажал ID)
"""
import sys

from Bio import SeqIO


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    inp, out = sys.argv[1], sys.argv[2]

    records = []
    for rec in SeqIO.parse(inp, "fasta"):
        for strand, nuc in ((+1, rec.seq), (-1, rec.seq.reverse_complement())):
            for frame in range(3):
                sub = nuc[frame:]
                sub = sub[: len(sub) - len(sub) % 3]
                prot = sub.translate(table=1, to_stop=False)
                rid = f"{rec.id}_frame{frame + 1}{'+' if strand == 1 else '-'}"
                records.append(f">{rid}\n{prot}")

    with open(out, "w") as fh:
        fh.write("\n".join(records) + "\n")


if __name__ == "__main__":
    main()
