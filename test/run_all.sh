#!/usr/bin/env bash
# =============================================================================
# run_all.sh — ПОЛНЫЙ цикл одной командой:
#   (setup, если ещё не делался) -> prepare_dataset -> run_gene для гена(ов) -> zip
#
# Использование:
#   ./run_all.sh <input.txt|.fasta> -g ГЕН          [опции]   # один ген
#   ./run_all.sh <input.txt|.fasta> --genes list.txt [опции]  # список генов
#   ./run_all.sh <input.txt|.fasta> --all            [опции]  # все гены файла
#
# Опции:
#   -n, --name ИМЯ   имя датасета (по умолчанию — имя файла без расширения)
#   --force          пересоздать датасет (делается и сам, если файл изменился)
#   --no-zip         не упаковывать результаты в архив
#   всё остальное (-e, -t, --min-pident, --top, ...) передаётся в run_gene.sh
#
# Результат: results/ДАТАСЕТ/ГЕН/..., results/ДАТАСЕТ/summary.tsv,
#            results/ДАТАСЕТ.zip
# =============================================================================
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

usage() { sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

INPUT=""; DATASET=""; FORCE=0; ZIP=1; ALL=0
GENES=(); GENES_FILE=""; PASS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -g|--gene)   GENES+=("${2:?нужен ген}"); shift 2 ;;
        --genes)     GENES_FILE="${2:?нужен файл}"; shift 2 ;;
        --all)       ALL=1; shift ;;
        -n|--name)   DATASET="${2:?нужно имя}"; shift 2 ;;
        --force)     FORCE=1; shift ;;
        --no-zip)    ZIP=0; shift ;;
        -h|--help)   usage; exit 0 ;;
        -*)          PASS+=("$1"); [[ $# -ge 2 && "$2" != -* ]] && { PASS+=("$2"); shift; }; shift ;;
        *)           [[ -z "$INPUT" ]] && INPUT="$1" || die "Лишний аргумент: $1"; shift ;;
    esac
done

[[ -n "$INPUT" ]] || { usage; exit 1; }
[[ -f "$INPUT" ]] || die "Файл не найден: $INPUT"
[[ ${#GENES[@]} -gt 0 || -n "$GENES_FILE" || $ALL -eq 1 ]] \
    || die "Укажите ген: -g ГЕН, --genes файл или --all"

# ---------- 0. Окружение (setup.sh идемпотентен) ----------
if [[ ! -f "$VENV_DIR/.deps_installed" ]]; then
    log "Окружение не настроено — запускаю setup.sh..."
    "$PROJECT_DIR/setup.sh" || die "setup.sh завершился с ошибкой"
fi

[[ -n "$DATASET" ]] || DATASET="$(dataset_name_from_file "$INPUT")"
DATASET="$(safe_name "$DATASET")"

# ---------- 1. Датасет (пересоздаём, если содержимое файла изменилось) ----------
HASH="$(sha256sum "$INPUT" | cut -d' ' -f1)"
STAMP="$PROJECT_DIR/data/processed/$DATASET/.source_sha256"
if [[ -f "$STAMP" && "$(cat "$STAMP")" != "$HASH" ]]; then
    log "Файл изменился с прошлой подготовки — пересоздаю датасет '$DATASET'."
    FORCE=1
fi

prep_args=("$INPUT" -n "$DATASET")
[[ $FORCE -eq 1 ]] && prep_args+=(--force)
"$PROJECT_DIR/prepare_dataset.sh" "${prep_args[@]}" || die "Подготовка датасета не удалась"
echo "$HASH" > "$STAMP"

CAND="$PROJECT_DIR/data/processed/$DATASET/candidates.fasta"

# ---------- 2. Список генов ----------
if [[ $ALL -eq 1 ]]; then
    mapfile -t GENES < <(grep '^>' "$CAND" | sed 's/^>//; s/[[:space:]].*//')
elif [[ -n "$GENES_FILE" ]]; then
    [[ -f "$GENES_FILE" ]] || die "Файл со списком генов не найден: $GENES_FILE"
    while IFS= read -r g || [[ -n "$g" ]]; do
        g="$(echo "${g//$'\r'/}" | xargs)"
        [[ -z "$g" || "$g" == \#* ]] && continue
        GENES+=("$g")
    done < "$GENES_FILE"
fi
[[ ${#GENES[@]} -gt 0 ]] || die "Список генов пуст"

# ---------- 3. Анализ ----------
RES="$PROJECT_DIR/results/$DATASET"
mkdir -p "$RES"
SUMMARY="$RES/summary.tsv"
printf 'gene\tstatus\tbest_frame\tn_hits\ttree\n' > "$SUMMARY"

n_ok=0; n_fail=0; i=0
for gene in "${GENES[@]}"; do
    i=$((i + 1))
    log "=== [$i/${#GENES[@]}] $gene ==="
    gs="$(safe_name "$gene")"
    d="$RES/$gs"
    rm -f "$d/${gs}_best_frame.txt" "$d/${gs}_blastp_hit_ids.txt"   # чтобы в сводку не попало старое
    if "$PROJECT_DIR/run_gene.sh" -d "$DATASET" -g "$gene" "${PASS[@]}"; then
        status=ok; n_ok=$((n_ok + 1))
    else
        status=failed; n_fail=$((n_fail + 1))
    fi
    bf="$(cat "$d/${gs}_best_frame.txt" 2>/dev/null || echo -)"
    nh="$( { wc -l < "$d/${gs}_blastp_hit_ids.txt"; } 2>/dev/null || echo 0)"
    tree="$([[ -s "$d/phylogeny/${gs}_blastp.treefile" ]] && echo yes || echo no)"
    printf '%s\t%s\t%s\t%s\t%s\n' "$gene" "$status" "$bf" "$nh" "$tree" >> "$SUMMARY"
done

# ---------- 4. Архив ----------
if [[ $ZIP -eq 1 ]]; then
    rm -f "$RES.zip"
    (cd "$PROJECT_DIR/results" && python3 -m zipfile -c "$DATASET.zip" "$DATASET") \
        && log "Архив: $RES.zip"
fi

log "Итог: успешно $n_ok, с ошибкой $n_fail. Сводка: $SUMMARY"
[[ $n_fail -eq 0 ]]
