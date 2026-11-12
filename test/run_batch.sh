#!/usr/bin/env bash
# =============================================================================
# run_batch.sh — запуск run_gene.sh для списка генов из файла.
#
# Использование:
#   ./run_batch.sh ДАТАСЕТ genes.txt [доп. опции для run_gene.sh]
#
# genes.txt: по одному идентификатору гена в строке (пустые строки и
#            строки, начинающиеся с #, игнорируются).
# Пример:
#   ./run_batch.sh supplementary_file_unannotated_alignments genes.txt -t 4 --top 30
# =============================================================================
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

[[ $# -ge 2 ]] || { sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1; }

DATASET="$1"; GENES_FILE="$2"; shift 2
[[ -f "$GENES_FILE" ]] || die "Файл со списком генов не найден: $GENES_FILE"

ok=(); failed=()
while IFS= read -r gene || [[ -n "$gene" ]]; do
    gene="${gene//$'\r'/}"
    gene="$(echo "$gene" | xargs)"
    [[ -z "$gene" || "$gene" == \#* ]] && continue

    log "=== $gene ==="
    if "$PROJECT_DIR/run_gene.sh" -d "$DATASET" -g "$gene" "$@"; then
        ok+=("$gene")
    else
        failed+=("$gene")
    fi
done < "$GENES_FILE"

log "Успешно: ${#ok[@]} | С ошибкой: ${#failed[@]}"
if [[ ${#failed[@]} -gt 0 ]]; then
    log "Не обработаны: ${failed[*]}"
    exit 1
fi
