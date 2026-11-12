#!/usr/bin/env python3
"""Извлекает из FASTA первую запись, чей ID содержит (или точно равен) запросу.

Использование:
  extract_gene.py --fasta in.fasta --query FBgn0262810 --out gene.fasta [--exact]

В stdout печатается полный ID найденной записи (его использует run_gene.sh).
Сообщения — в stderr.
"""
import argparse
import sys

from Bio import SeqIO


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--query", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--exact", action="store_true", help="точное совпадение ID")
    args = ap.parse_args()

    exact, partial = [], []
    for rec in SeqIO.parse(args.fasta, "fasta"):
        if rec.id == args.query:
            exact.append(rec)
        elif not args.exact and args.query in rec.id:
            partial.append(rec)
    # Точное совпадение ID всегда приоритетнее подстроки
    # (иначе запрос "gene1" мог бы взять "gene10", если тот стоит раньше в файле)
    matches = exact or partial

    if not matches:
        sys.exit(f"ОШИБКА: запись '{args.query}' не найдена в {args.fasta}")

    rec = matches[0]
    if len(matches) > 1:
        print(
            f"ВНИМАНИЕ: найдено {len(matches)} записей, использую первую: {rec.id}",
            file=sys.stderr,
        )

    SeqIO.write(rec, args.out, "fasta")
    print(f"Извлечён: {rec.id}, длина: {len(rec.seq)}", file=sys.stderr)
    print(rec.id)


if __name__ == "__main__":
    main()
