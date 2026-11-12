#!/usr/bin/env bash
# =============================================================================
# run_gene.sh — анализ ОДНОГО гена (шаги 7, 9, 11–22):
#   извлечение -> blastn -> 6 рамок -> blastp -> фильтрация -> MSA -> дерево
# Можно запускать сколько угодно раз для разных генов и датасетов.
#
# Использование:
#   ./run_gene.sh -d ДАТАСЕТ -g ГЕН [опции]
#
# Обязательные:
#   -d, --dataset   имя датасета (см. ./prepare_dataset.sh)
#   -g, --gene      идентификатор гена (подстрока заголовка, напр. FBgn0262810)
# Опции (в скобках — значения по умолчанию):
#   -e, --evalue        порог e-value               (1e-5)
#   -t, --threads       число потоков               (2)
#   --min-pident        мин. % идентичности         (30)
#   --min-qcov          мин. покрытие запроса, %    (50)
#   --min-len           мин. длина выравнивания     (30)
#   --top               сколько лучших хитов брать  (50)
#   --max-targets       max_target_seqs для BLAST   (100)
#   --min-orf           мин. длина ORF для getorf   (90)
#   --orf-find          режим -find для getorf      (3)
#   -h, --help          справка
#
# Результаты: results/ДАТАСЕТ/ГЕН/   (логи: logs/)
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

usage() {
    sed -n '2,27p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

DATASET=""; GENE=""
EVALUE="1e-5"; THREADS=2
MIN_PIDENT=30; MIN_QCOV=50; MIN_LEN=30
TOP_N=50; MAX_TARGETS=100
MIN_ORF=90; ORF_FIND=3

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dataset)    DATASET="${2:?}"; shift 2 ;;
        -g|--gene)       GENE="${2:?}"; shift 2 ;;
        -e|--evalue)     EVALUE="${2:?}"; shift 2 ;;
        -t|--threads)    THREADS="${2:?}"; shift 2 ;;
        --min-pident)    MIN_PIDENT="${2:?}"; shift 2 ;;
        --min-qcov)      MIN_QCOV="${2:?}"; shift 2 ;;
        --min-len)       MIN_LEN="${2:?}"; shift 2 ;;
        --top)           TOP_N="${2:?}"; shift 2 ;;
        --max-targets)   MAX_TARGETS="${2:?}"; shift 2 ;;
        --min-orf)       MIN_ORF="${2:?}"; shift 2 ;;
        --orf-find)      ORF_FIND="${2:?}"; shift 2 ;;
        -h|--help)       usage; exit 0 ;;
        *)               die "Неизвестный аргумент: $1 (см. --help)" ;;
    esac
done

[[ -n "$DATASET" && -n "$GENE" ]] || { usage; exit 1; }

activate_venv
require_cmd blastn blastp blastdbcmd mafft getorf
IQTREE="$(find_iqtree)" || die "IQ-TREE не найден. Выполните ./setup.sh"

# ---------- Пути ----------
GENE_SAFE="$(safe_name "$GENE")"
CAND="$PROJECT_DIR/data/processed/$DATASET/candidates.fasta"
BLASTN_DB="$PROJECT_DIR/databases/blastn/$DATASET/cds"
BLASTP_DB="$PROJECT_DIR/databases/blastp/$DATASET/proteins"

if [[ ! -s "$CAND" ]]; then
    avail="$(ls "$PROJECT_DIR/data/processed" 2>/dev/null | tr '\n' ' ')"
    die "Датасет '$DATASET' не подготовлен. Доступные: ${avail:-нет}. Запустите ./prepare_dataset.sh"
fi

OUT="$PROJECT_DIR/results/$DATASET/$GENE_SAFE"
PHYLO="$OUT/phylogeny"
mkdir -p "$OUT" "$PHYLO" "$PROJECT_DIR/logs"

LOG="$PROJECT_DIR/logs/${DATASET}_${GENE_SAFE}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

log "Датасет: $DATASET | Ген: $GENE | Потоки: $THREADS | e-value: $EVALUE"

P="$OUT/$GENE_SAFE"   # префикс выходных файлов

# ---------- 7. Извлечение последовательности гена ----------
log "[1/10] Извлекаю последовательность гена..."
REC_ID="$(python3 "$SCRIPTS_DIR/extract_gene.py" \
            --fasta "$CAND" --query "$GENE" --out "$P.fasta")"

# ---------- 9. blastn ----------
log "[2/10] blastn..."
blastn \
    -query "$P.fasta" \
    -db "$BLASTN_DB" \
    -evalue "$EVALUE" \
    -max_target_seqs "$MAX_TARGETS" \
    -max_hsps 1 \
    -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore stitle" \
    -num_threads "$THREADS" \
    -out "${P}_blastn.tsv"

# ---------- 11. Шесть рамок считывания ----------
log "[3/10] Трансляция в 6 рамках..."
python3 "$SCRIPTS_DIR/six_frames.py" "$P.fasta" "${P}_6frames.faa"

# ---------- 12. blastp ----------
log "[4/10] blastp..."
blastp \
    -query "${P}_6frames.faa" \
    -db "$BLASTP_DB" \
    -evalue "$EVALUE" \
    -max_target_seqs "$MAX_TARGETS" \
    -max_hsps 1 \
    -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore stitle" \
    -num_threads "$THREADS" \
    -out "${P}_blastp.tsv"

[[ -s "${P}_blastp.tsv" ]] || die "blastp не вернул результатов для $GENE"

# ---------- 13. Фильтрация (self-хиты — белки самого гена — исключаются) ----------
log "[5/10] Фильтрую хиты..."
awk -F'\t' -v OFS='\t' \
    -v pid="$MIN_PIDENT" -v qcov="$MIN_QCOV" -v ev="$EVALUE" -v len="$MIN_LEN" \
    -v self="${REC_ID}_" \
    '$3+0 >= pid && $5+0 >= qcov && $6+0 <= ev+0 && $4+0 >= len \
     && substr($2, 1, length(self)) != self' \
    "${P}_blastp.tsv" > "${P}_blastp.filtered.tsv"

[[ -s "${P}_blastp.filtered.tsv" ]] \
    || die "После фильтрации не осталось хитов. Ослабьте --min-pident / --min-qcov / --evalue"

# ---------- 17 (автоматически). Лучшая рамка = запрос с макс. bitscore ----------
# Считаем по отфильтрованным хитам БЕЗ self-хитов: иначе каждая рамка
# находит саму себя со 100% идентичности, и «побеждает» просто самая длинная
# рамка без стоп-кодонов, а не та, что лучше всего похожа на другие гены.
BEST_FRAME="$(awk -F'\t' '$7+0>m {m=$7+0; q=$1} END{print q}' "${P}_blastp.filtered.tsv")"
log "      Лучшая рамка: $BEST_FRAME"
echo "$BEST_FRAME" > "${P}_best_frame.txt"

# ---------- 14. Топ-N по bitscore ----------
sort -t$'\t' -k7,7gr "${P}_blastp.filtered.tsv" -o "${P}_blastp.sorted.tmp"
head -n "$TOP_N" "${P}_blastp.sorted.tmp" > "${P}_blastp.top${TOP_N}.tsv"
rm -f "${P}_blastp.sorted.tmp"

# ---------- 15. Уникальные ID хитов ----------
cut -f2 "${P}_blastp.top${TOP_N}.tsv" | sort -u > "${P}_blastp_hit_ids.txt"
log "      Уникальных хитов: $(wc -l < "${P}_blastp_hit_ids.txt")"

# ---------- 16. Последовательности хитов ----------
log "[6/10] Извлекаю белки-хиты..."
blastdbcmd -db "$BLASTP_DB" -entry_batch "${P}_blastp_hit_ids.txt" -out "${P}_blastp_hits.faa"

# ---------- 17. Белок в лучшей рамке ----------
python3 "$SCRIPTS_DIR/extract_gene.py" \
    --fasta "${P}_6frames.faa" --query "$BEST_FRAME" --exact \
    --out "${P}_best_frame.faa" >/dev/null

# ---------- 18. Вход для выравнивания ----------
cat "${P}_best_frame.faa" "${P}_blastp_hits.faa" > "${P}_msa_input.faa"

# ---------- 19. MAFFT ----------
log "[7/10] Множественное выравнивание (MAFFT)..."
mafft --auto --thread "$THREADS" "${P}_msa_input.faa" > "${P}_alignment.fasta"

# ---------- 20. IQ-TREE ----------
N_SEQ="$(count_fasta "${P}_alignment.fasta")"
if [[ "$N_SEQ" -ge 4 ]]; then
    log "[8/10] Филогения (IQ-TREE, $N_SEQ последовательностей)..."
    "$IQTREE" \
        -s "${P}_alignment.fasta" \
        -m MFP \
        -B 1000 \
        --alrt 1000 \
        -T "$THREADS" \
        --prefix "$PHYLO/${GENE_SAFE}_blastp" \
        --redo
else
    log "[8/10] Пропускаю IQ-TREE: нужно ≥ 4 последовательностей, найдено $N_SEQ."
fi

# ---------- 21. ORF-ы ----------
log "[9/10] Поиск ORF (getorf)..."
getorf \
    -sequence "$P.fasta" \
    -outseq "${P}_orfs.fasta" \
    -minsize "$MIN_ORF" \
    -find "$ORF_FIND"

# ---------- 22. Отрисовка дерева ----------
TREEFILE="$PHYLO/${GENE_SAFE}_blastp.treefile"
if [[ -s "$TREEFILE" ]]; then
    log "[10/10] Рисую дерево..."
    python3 "$SCRIPTS_DIR/draw_tree.py" "$TREEFILE" "$PHYLO/${GENE_SAFE}_blastp"
else
    log "[10/10] Дерево не построено — пропускаю отрисовку."
fi

log "Готово: $OUT"
log "Лог: $LOG"
