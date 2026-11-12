#!/usr/bin/env python3
"""Конвертирует текстовый файл с блоками 'raw sequence for <ID>' в FASTA.

Использование: txt_to_fasta.py <input.txt> <output.fasta>
"""
import re
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    inp, out = sys.argv[1], sys.argv[2]

    with open(inp, encoding="utf-8") as fh:
        text = fh.read()

    pattern = re.compile(
        r"raw sequence for\s+(\S+)(.*?)(?=raw sequence for\s+\S+|\Z)",
        re.IGNORECASE | re.DOTALL,
    )

    records = []
    for m in pattern.finditer(text):
        rid, body = m.group(1), m.group(2)

        # Если внутри блока есть строка заголовка FASTA (">") — берём всё после неё
        header = re.search(r"^>.*$", body, re.MULTILINE)
        seq_part = body[header.end():] if header else body

        # Оставляем только строки, состоящие из ACGTN (и пробелов)
        seq_lines = []
        for line in seq_part.splitlines():
            line = line.strip()
            if line and re.fullmatch(r"[ACGTNacgtn\s]+", line):
                seq_lines.append(re.sub(r"\s+", "", line))

        seq = "".join(seq_lines).upper()
        if seq:
            records.append(SeqRecord(Seq(seq), id=rid, description=""))

    SeqIO.write(records, out, "fasta")
    print(f"Записано последовательностей: {len(records)}", file=sys.stderr)


if __name__ == "__main__":
    main()
