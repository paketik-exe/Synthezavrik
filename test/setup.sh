#!/usr/bin/env bash
# =============================================================================
# setup.sh — РАЗОВАЯ установка (шаги 1–3 исходной инструкции).
# Безопасно запускать повторно: уже установленное пропускается.
#
# Использование:  ./setup.sh [--force]
#   --force  переустановить Python-пакеты в виртуальном окружении
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

FORCE=0
[[ "${1:-}" == "--force" ]] && FORCE=1

# ---------- 1. Системные пакеты ----------
APT_PACKAGES=(ncbi-blast+ mafft emboss iqtree python3 python3-venv python3-pip)

missing=0
for c in blastn blastp makeblastdb blastdbcmd mafft transeq getorf python3; do
    command -v "$c" >/dev/null 2>&1 || { log "Не найдено: $c"; missing=1; }
done
find_iqtree >/dev/null 2>&1 || { log "Не найдено: iqtree"; missing=1; }
python3 -c "import venv, ensurepip" >/dev/null 2>&1 || { log "Не найден python3-venv"; missing=1; }

if [[ $missing -eq 1 ]]; then
    log "Устанавливаю системные пакеты..."
    SUDO=""
    [[ $EUID -ne 0 ]] && SUDO="sudo"
    $SUDO apt-get update
    DEBIAN_FRONTEND=noninteractive $SUDO apt-get install -y "${APT_PACKAGES[@]}"
else
    log "Системные пакеты уже установлены — пропускаю."
fi

IQTREE_BIN="$(find_iqtree)" || die "IQ-TREE не установился"
if ! "$IQTREE_BIN" --version 2>&1 | grep -qi 'version 2'; then
    log "ВНИМАНИЕ: пайплайн рассчитан на IQ-TREE 2.x (ключи --alrt, --prefix)."
    log "          Установлена другая версия — установите IQ-TREE 2 (conda/bioconda или с сайта)."
fi

# ---------- 2. Виртуальное окружение Python ----------
if [[ ! -f "$VENV_DIR/bin/activate" ]]; then
    log "Создаю виртуальное окружение..."
    python3 -m venv "$VENV_DIR"
    FORCE=1
fi

if [[ $FORCE -eq 1 || ! -f "$VENV_DIR/.deps_installed" ]]; then
    log "Устанавливаю Python-пакеты..."
    # shellcheck disable=SC1091
    source "$VENV_DIR/bin/activate"
    pip install --upgrade pip
    pip install biopython pandas numpy scipy pyyaml streamlit plotly \
                matplotlib jinja2 pydantic pytest
    touch "$VENV_DIR/.deps_installed"
else
    log "Python-пакеты уже установлены — пропускаю."
fi

# ---------- 3. Структура каталогов ----------
mkdir -p "$PROJECT_DIR"/data/{uploads,processed,reference/{genomes,cds,proteins,annotations,phylogeny},test}
mkdir -p "$PROJECT_DIR"/databases/{blastn,blastp,metadata}
mkdir -p "$PROJECT_DIR"/results "$PROJECT_DIR"/scripts "$PROJECT_DIR"/logs

log "Готово. Дальше:"
log "  ./prepare_dataset.sh <файл_с_последовательностями.txt>"
log "  ./run_gene.sh -d <датасет> -g <ген>"
