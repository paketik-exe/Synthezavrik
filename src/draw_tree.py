#!/usr/bin/env python3
"""Рисует дерево: ASCII в консоль и файл, плюс PNG (если доступен matplotlib).

Использование: draw_tree.py <tree.treefile> <префикс_выходных_файлов>
Создаёт: <префикс>.ascii.txt и <префикс>.png
"""
import io
import sys
from contextlib import redirect_stdout

from Bio import Phylo


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    treefile, prefix = sys.argv[1], sys.argv[2]

    tree = Phylo.read(treefile, "newick")

    buf = io.StringIO()
    with redirect_stdout(buf):
        Phylo.draw_ascii(tree)
    ascii_tree = buf.getvalue()
    print(ascii_tree)
    with open(f"{prefix}.ascii.txt", "w") as fh:
        fh.write(ascii_tree)

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        n_leaves = tree.count_terminals()
        fig = plt.figure(figsize=(10, max(6, n_leaves * 0.25)))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=ax, do_show=False)
        fig.savefig(f"{prefix}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"PNG сохранён: {prefix}.png", file=sys.stderr)
    except Exception as exc:  # отрисовка PNG — необязательная часть
        print(f"PNG не создан: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()
