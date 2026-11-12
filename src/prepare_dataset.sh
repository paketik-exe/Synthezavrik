#!/usr/bin/env bash
# =============================================================================
# prepare_dataset.sh — подготовка ОДНОГО входного файла (шаги 5, 6, 8, 10):
#   txt -> candidates.fasta -> BLAST-базы (нуклеотидная и белковая)
# Запускается один раз на файл; повторно — только с --force.
#
# Использование:
#   ./prepare_dataset.sh <input.txt> [-n имя_датасета] [--force]
#
# Результат (для датасета DS):
#   data/processed/DS/candidates.fasta
#   data/processed/DS/candidates_proteins.fasta
#   databases/blastn/DS/cds.*      databases/blastp/DS/proteins.*
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

usage() {
    sed -n '2,13p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

INPUT=""
DATASET=""
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--name)  DATASET="${2:?нужно имя}"; shift 2 ;;
        --force)    FORCE=1; shift ;;
        -h|--help)  usage; exit 0 ;;
        -*)         die "Неизвестный аргумент: $1" ;;
        *)          [[ -z "$INPUT" ]] && INPUT="$1" || die "Лишний аргумент: $1"; shift ;;
    esac
done

[[ -n "$INPUT" ]] || { usage; exit 1; }
[[ -f "$INPUT" ]] || die "Файл не найден: $INPUT"

activate_venv
require_cmd makeblastdb transeq

[[ -n "$DATASET" ]] || DATASET="$(dataset_name_from_file "$INPUT")"
DATASET="$(safe_name "$DATASET")"

PROC_DIR="$PROJECT_DIR/data/processed/$DATASET"
CAND="$PROC_DIR/candidates.fasta"
PROT="$PROC_DIR/candidates_proteins.fasta"
BLASTN_DB="$PROJECT_DIR/databases/blastn/$DATASET/cds"
BLASTP_DB="$PROJECT_DIR/databases/blastp/$DATASET/proteins"

if [[ -s "$CAND" && -s "${BLASTN_DB}.nhr" && -s "${BLASTP_DB}.phr" && $FORCE -eq 0 ]]; then
    log "Датасет '$DATASET' уже подготовлен — пропускаю (используйте --force для пересоздания)."
    exit 0
fi

mkdir -p "$PROC_DIR" "$(dirname "$BLASTN_DB")" "$(dirname "$BLASTP_DB")" \
         "$PROJECT_DIR/data/uploads" "$PROJECT_DIR/logs"

# ---------- 5. Копируем исходник в uploads (если он ещё не там) ----------
UPLOAD="$PROJECT_DIR/data/uploads/$(basename "$INPUT")"
if [[ "$(realpath "$INPUT")" != "$(realpath -m "$UPLOAD")" ]]; then
    cp -f "$INPUT" "$UPLOAD"
    log "Файл скопирован в $UPLOAD"
fi

# ---------- 6. txt -> candidates.fasta ----------
log "Конвертирую в FASTA..."
python3 "$SCRIPTS_DIR/txt_to_fasta.py" "$UPLOAD" "$CAND"
[[ -s "$CAND" ]] || die "Не удалось извлечь последовательности из $INPUT"

# ---------- 8. Нуклеотидная BLAST-база ----------
log "Строю нуклеотидную BLAST-базу..."
makeblastdb -in "$CAND" -dbtype nucl -parse_seqids -out "$BLASTN_DB"

# ---------- 10. Трансляция в 6 рамках + белковая BLAST-база ----------
log "Транслирую в 6 рамках (transeq)..."
transeq -sequence "$CAND" -outseq "$PROT" -frame 6

log "Строю белковую BLAST-базу..."
makeblastdb -in "$PROT" -dbtype prot -parse_seqids -out "$BLASTP_DB"

log "Датасет '$DATASET' готов: $(count_fasta "$CAND") последовательностей."
log "Запуск анализа:  ./run_gene.sh -d $DATASET -g <ID_гена>"
