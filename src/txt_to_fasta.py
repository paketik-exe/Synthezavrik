#!/usr/bin/env python3
"""Конвертирует текстовый файл с блоками 'raw sequence for <ID>' в FASTA.

Если таких блоков нет, а файл уже в формате FASTA, — просто нормализует его.

Использование: txt_to_fasta.py <input.txt> <output.fasta>
"""
import re
import sys

BLOCK_RE = re.compile(
    r"raw sequence for\s+(\S+)(.*?)(?=raw sequence for\s+\S+|\Z)",
    re.IGNORECASE | re.DOTALL,
)
SEQ_LINE_RE = re.compile(r"[ACGTNacgtn\s]+")


def _seq_from_body(body: str) -> str:
    # Если внутри блока есть строка заголовка FASTA (">") — берём всё после неё
    header = re.search(r"^>.*$", body, re.MULTILINE)
    seq_part = body[header.end():] if header else body

    # Оставляем только строки, состоящие из ACGTN (и пробелов)
    seq_lines = []
    for line in seq_part.splitlines():
        line = line.strip()
        if line and SEQ_LINE_RE.fullmatch(line):
            seq_lines.append(re.sub(r"\s+", "", line))
    return "".join(seq_lines).upper()


def _parse_fasta(text: str) -> list[tuple[str, str]]:
    records, rid, buf = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if rid is not None:
                records.append((rid, "".join(buf).upper()))
            parts = line[1:].split()
            rid, buf = (parts[0] if parts else f"seq{len(records) + 1}"), []
        elif rid is not None and line:
            buf.append(re.sub(r"\s+", "", line))
    if rid is not None:
        records.append((rid, "".join(buf).upper()))
    return records


def _dedupe(records: list[tuple[str, str]]) -> list[tuple[str, str]]:
    """makeblastdb -parse_seqids падает на повторяющихся ID — переименовываем."""
    seen: dict[str, int] = {}
    out = []
    for rid, seq in records:
        if rid in seen:
            seen[rid] += 1
            new_id = f"{rid}_dup{seen[rid]}"
            print(f"ВНИМАНИЕ: повтор ID '{rid}' -> '{new_id}'", file=sys.stderr)
            rid = new_id
        else:
            seen[rid] = 1
        out.append((rid, seq))
    return out


def parse_records(text: str) -> list[tuple[str, str]]:
    """Возвращает [(id, последовательность)] из текста в любом из двух форматов."""
    records = [(m.group(1), _seq_from_body(m.group(2))) for m in BLOCK_RE.finditer(text)]
    if not records and text.lstrip().startswith(">"):
        records = _parse_fasta(text)
    records = [(rid, seq) for rid, seq in records if seq]
    return _dedupe(records)


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    inp, out = sys.argv[1], sys.argv[2]

    with open(inp, encoding="utf-8", errors="replace") as fh:
        records = parse_records(fh.read())

    with open(out, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    print(f"Записано последовательностей: {len(records)}", file=sys.stderr)


if __name__ == "__main__":
    main()
