#!/usr/bin/env bash
# run_app.sh — запуск веб-интерфейса (http://localhost:8501).
# При первом запуске сам выполнит setup.sh.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

[[ -f "$VENV_DIR/.deps_installed" ]] || "$PROJECT_DIR/setup.sh"
activate_venv
exec streamlit run "$PROJECT_DIR/app.py" --server.maxUploadSize 1024 "$@"
