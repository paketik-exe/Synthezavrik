#!/usr/bin/env python3
"""Выбирает самый длинный ORF гена по выводу getorf.

Использование:
  pick_orf.py <gene.fasta> <orfs_aa.fasta (getorf -find 1)> <orfs_nt.fasta (getorf -find 3)> <outprefix>
              [--protein known.faa]

Создаёт:
  <outprefix>.faa  — белок ORF (запрос для blastp/tblastn)
  <outprefix>.fna  — нуклеотиды ORF (запрос для blastn)
  <outprefix>.tsv  — метрики: orf_id, orf_len_aa, orf_len_nt, strand, start, end, gene_len, orf_frac, n_orfs
С --protein берётся известный белок, а для blastn — весь ген.
"""
import re
import sys


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def main() -> None:
    args = sys.argv[1:]
    protein = None
    if "--protein" in args:
        i = args.index("--protein")
        protein = args[i + 1]
        del args[i:i + 2]
    if len(args) != 4:
        sys.exit(__doc__)
    gene_fa, orf_aa, orf_nt, prefix = args

    gene_name, gene_seq = read_fasta(gene_fa)[0]
    gene_len = len(gene_seq)
    aa = read_fasta(orf_aa)
    nt = {h.split()[0]: s for h, s in read_fasta(orf_nt)}

    if protein:
        pname, pseq = read_fasta(protein)[0]
        pseq = pseq.rstrip("*")
        orf_id, aa_seq, nt_seq, strand, start, end = pname.split()[0], pseq, gene_seq, "?", 1, gene_len
    else:
        if not aa:
            aa_seq = ""
            orf_id, nt_seq, strand, start, end = "-", gene_seq, "?", 0, 0
        else:
            head, aa_seq = max(aa, key=lambda r: len(r[1]))
            orf_id = head.split()[0]
            m = re.search(r"\[(\d+)\s*-\s*(\d+)\]", head)
            start, end = (int(m.group(1)), int(m.group(2))) if m else (0, 0)
            strand = "-" if "REVERSE" in head.upper() else "+"
            nt_seq = nt.get(orf_id, gene_seq)
            aa_seq = aa_seq.rstrip("*")

    with open(prefix + ".faa", "w") as fh:
        fh.write(f">{orf_id}\n{aa_seq}\n")
    with open(prefix + ".fna", "w") as fh:
        fh.write(f">{orf_id}\n{nt_seq}\n")
    frac = round(100.0 * len(aa_seq) * 3 / gene_len, 1) if gene_len else 0
    cols = ["orf_id", "orf_len_aa", "orf_len_nt", "strand", "start", "end", "gene_len", "orf_frac", "n_orfs"]
    vals = [orf_id, len(aa_seq), len(aa_seq) * 3, strand, start, end, gene_len, frac, len(aa)]
    with open(prefix + ".tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n" + "\t".join(map(str, vals)) + "\n")
    print(f"ORF: {orf_id}, {len(aa_seq)} а.к. (цепь {strand}), "
          f"{frac}% длины гена; всего ORF: {len(aa)}", file=sys.stderr)


if __name__ == "__main__":
    main()
