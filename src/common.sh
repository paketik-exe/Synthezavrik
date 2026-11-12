#!/usr/bin/env bash
# Общие функции и пути. Подключается через: source common.sh
# Не запускается напрямую.

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$PROJECT_DIR/.venv"
# Python-скрипты могут лежать в scripts/ или рядом с .sh — берём, где они есть
if [[ -f "$PROJECT_DIR/scripts/txt_to_fasta.py" ]]; then
    SCRIPTS_DIR="$PROJECT_DIR/scripts"
else
    SCRIPTS_DIR="$PROJECT_DIR"
fi

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }
die() { log "ОШИБКА: $*"; exit 1; }

activate_venv() {
    [[ -f "$VENV_DIR/bin/activate" ]] \
        || die "Окружение не найдено. Сначала выполните: ./setup.sh"
    # shellcheck disable=SC1091
    source "$VENV_DIR/bin/activate"
}

require_cmd() {
    for c in "$@"; do
        command -v "$c" >/dev/null 2>&1 || die "Не найдена программа '$c'. Выполните ./setup.sh"
    done
}

# iqtree в разных дистрибутивах называется по-разному
find_iqtree() {
    if command -v iqtree2 >/dev/null 2>&1; then echo iqtree2
    elif command -v iqtree >/dev/null 2>&1; then echo iqtree
    else return 1
    fi
}

# Имя датасета по имени входного файла: /path/my.file.txt -> my.file
dataset_name_from_file() {
    local base
    base="$(basename "$1")"
    echo "${base%.*}"
}

# Безопасное для имён файлов представление идентификатора
safe_name() { echo "${1//[^A-Za-z0-9._-]/_}"; }

count_fasta() { grep -c '^>' "$1" || true; }
