#!/usr/bin/env bash
# =============================================================================
# classify_gene.sh — эволюционное происхождение гена (алгоритм «Синтезаврики»):
#   1. Есть гомологи?  нет → 3. синтенное некодирующее выравнивание у аутгруппы → De novo
#                                 иначе связь с TE/ncRNA → Иное происхождение, иначе Неопределённо
#                      да  → 2.1 в том же геноме → Дупликация
#                            2.2 у далёких видов в синтенной позиции → Древний ген
#                            2.3 вне синтении: мозаика → Слияние/перестройка;
#                                другая ветвь жизни → HGT; иначе Неопределённо
#
# Использование:
#   ./classify_gene.sh -d ДАТАСЕТ -g ГЕН [-g ГЕН2 ...] [опции]
#   ./classify_gene.sh -d ДАТАСЕТ --all [опции]
#
# Виды:
#   --focal ВИД            вид гена (по умолчанию — префикс ID до "_")
#   --ingroup a,b          близкие виды (ортологи у них не мешают de novo)
#   --outgroup a,b         далёкие виды / аутгруппа (по умолчанию — все прочие)
#   --species-regex RE     вид из ID (по умолчанию '^([^_]+)_')
#   --gene-key-regex RE    общий ключ гена у разных видов (по умолчанию '(FBgn\d+)');
#                          записи датасета с тем же ключом считаются синтенными
# Данные (все необязательны, но каждая включает свою ветку алгоритма):
#   --ref-proteins DIR     протеомы видов: DIR/<вид>.faa            (2.1, 2.2)
#   --ref-genomes DIR      геномы видов: DIR/<вид>.fa               (2.2, 3, синтения)
#   --focal-genome FILE    геном фокального вида                     (синтения)
#   --focal-gff FILE       его аннотация GFF3 (координаты гена)      (синтения)
#   --locus CHR:START-END  координаты гена вместо GFF (для одного гена)
#   --flank N              длина фланков для синтении, нт (10000)
#   --hgt-db DIR           белки другой ветви жизни: DIR/<группа>.faa  (2.3.2)
#   --te-db FILE           библиотека мобильных элементов (FASTA)    (3.2)
#   --ncrna-db FILE        последовательности ncRNA (FASTA)          (3.2)
#   --protein FILE         известный белок гена (иначе — самый длинный ORF)
#   --online               искать онлайн в NCBI, ничего не скачивая (см. ниже)
#   --online-scope ENTREZ  где искать белковые гомологи онлайн ('"Insecta"[Organism]')
#   --genome-species a,b   геномы для онлайн-поиска некодирующих участков и синтении
#                          (по умолчанию = --outgroup, иначе dyak,dere,dana,dpse,dvir)
#   --refresh              не брать результаты онлайн-поиска из кэша
#   --use-references       взять всё, что скачал ./fetch_references.sh (data/reference/)
#   --tree                 построить дерево (run_gene.sh), если его нет
#   --require-synteny      непроверенная синтения считается отсутствующей
# Пороги (по умолчанию):
#   -e 1e-5  --min-pident 30  --min-cov 50  --part-cov 20  --max-overlap 20
#   --min-orf-aa 30  --nt-pident 70  --nt-cov 50  --min-ai 45
#   --te-evalue 1e-10  --te-pident 80  --te-cov 50  --min-uf 95  --min-alrt 80
#   --max-targets 500  -t 2
#
# Онлайн-режим (--online): удалённый BLAST на серверах NCBI + API NCBI Gene/Orthologs:
#   blastp  → refseq_protein (гомологи в --online-scope; другие ветви жизни для HGT)
#   tblastn/blastn → refseq_genomic (геномы --genome-species: некодирующие гомологи)
#   blastn фланков гена → refseq_genomic (синтенные окна); ортологи NCBI (синтения)
#   blastn гена → refseq_rna, только ncRNA; мобильные элементы — набор Bergman lab (~1 МБ)
#   Каждый удалённый поиск идёт от минуты до 10+ минут; результаты кэшируются.
#   В онлайн-режиме «далёкие виды» = все виды, кроме фокального и --ingroup.
#
# Результаты: results/ДАТАСЕТ/ГЕН/origin/{origin.tsv, origin_report.md, *.tsv}
#             results/ДАТАСЕТ/origin_summary.tsv
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

usage() { sed -n '2,60p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

DATASET=""; GENES=(); ALL=0
FOCAL=""; INGROUP=""; OUTGROUP=""; SPECIES_RE='^([^_]+)_'; KEY_RE='(FBgn\d+)'
REF_PROT=""; REF_GEN=""; FOCAL_GENOME=""; FOCAL_GFF=""; LOCUS=""; FLANK=10000
HGT_DB=""; TE_DB=""; NCRNA_DB=""; PROTEIN=""; BUILD_TREE=0; REQ_SYN=0; USE_REFS=0
ONLINE=0; REFRESH=0; ONLINE_SCOPE='"Insecta"[Organism]'; GEN_SPECIES=""
EVALUE="1e-5"; MIN_PIDENT=30; MIN_COV=50; PART_COV=20; MAX_OVERLAP=20
MIN_ORF_AA=30; NT_PIDENT=70; NT_COV=50; MIN_AI=45
TE_EVALUE="1e-10"; TE_PIDENT=80; TE_COV=50; MIN_UF=95; MIN_ALRT=80
MAX_TARGETS=500; THREADS=2

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dataset)      DATASET="${2:?}"; shift 2 ;;
        -g|--gene)         GENES+=("${2:?}"); shift 2 ;;
        --all)             ALL=1; shift ;;
        --focal)           FOCAL="${2:?}"; shift 2 ;;
        --ingroup)         INGROUP="${2:?}"; shift 2 ;;
        --outgroup)        OUTGROUP="${2:?}"; shift 2 ;;
        --species-regex)   SPECIES_RE="${2:?}"; shift 2 ;;
        --gene-key-regex)  KEY_RE="${2:?}"; shift 2 ;;
        --ref-proteins)    REF_PROT="${2:?}"; shift 2 ;;
        --ref-genomes)     REF_GEN="${2:?}"; shift 2 ;;
        --focal-genome)    FOCAL_GENOME="${2:?}"; shift 2 ;;
        --focal-gff)       FOCAL_GFF="${2:?}"; shift 2 ;;
        --locus)           LOCUS="${2:?}"; shift 2 ;;
        --flank)           FLANK="${2:?}"; shift 2 ;;
        --hgt-db)          HGT_DB="${2:?}"; shift 2 ;;
        --te-db)           TE_DB="${2:?}"; shift 2 ;;
        --ncrna-db)        NCRNA_DB="${2:?}"; shift 2 ;;
        --protein)         PROTEIN="${2:?}"; shift 2 ;;
        --tree)            BUILD_TREE=1; shift ;;
        --use-references)  USE_REFS=1; shift ;;
        --online)          ONLINE=1; shift ;;
        --online-scope)    ONLINE_SCOPE="${2:?}"; shift 2 ;;
        --genome-species)  GEN_SPECIES="${2:?}"; shift 2 ;;
        --refresh)         REFRESH=1; shift ;;
        --require-synteny) REQ_SYN=1; shift ;;
        -e|--evalue)       EVALUE="${2:?}"; shift 2 ;;
        --min-pident)      MIN_PIDENT="${2:?}"; shift 2 ;;
        --min-cov)         MIN_COV="${2:?}"; shift 2 ;;
        --part-cov)        PART_COV="${2:?}"; shift 2 ;;
        --max-overlap)     MAX_OVERLAP="${2:?}"; shift 2 ;;
        --min-orf-aa)      MIN_ORF_AA="${2:?}"; shift 2 ;;
        --nt-pident)       NT_PIDENT="${2:?}"; shift 2 ;;
        --nt-cov)          NT_COV="${2:?}"; shift 2 ;;
        --min-ai)          MIN_AI="${2:?}"; shift 2 ;;
        --te-evalue)       TE_EVALUE="${2:?}"; shift 2 ;;
        --te-pident)       TE_PIDENT="${2:?}"; shift 2 ;;
        --te-cov)          TE_COV="${2:?}"; shift 2 ;;
        --min-uf)          MIN_UF="${2:?}"; shift 2 ;;
        --min-alrt)        MIN_ALRT="${2:?}"; shift 2 ;;
        --max-targets)     MAX_TARGETS="${2:?}"; shift 2 ;;
        -t|--threads)      THREADS="${2:?}"; shift 2 ;;
        -h|--help)         usage; exit 0 ;;
        *)                 die "Неизвестный аргумент: $1 (см. --help)" ;;
    esac
done

[[ -n "$DATASET" && ( ${#GENES[@]} -gt 0 || $ALL -eq 1 ) ]] || { usage; exit 1; }

# ---------- Референсы из fetch_references.sh (явные ключи важнее) ----------
if [[ $USE_REFS -eq 1 ]]; then
    R="$PROJECT_DIR/data/reference"
    has_fa() { [[ -d "$1" ]] && compgen -G "$1/*.fa*" >/dev/null; }
    [[ -n "$REF_PROT" ]]     || { has_fa "$R/proteins" && REF_PROT="$R/proteins"; } || true
    [[ -n "$REF_GEN" ]]      || { has_fa "$R/genomes"  && REF_GEN="$R/genomes"; } || true
    [[ -n "$HGT_DB" ]]       || { has_fa "$R/hgt"      && HGT_DB="$R/hgt"; } || true
    [[ -n "$FOCAL_GENOME" ]] || { [[ -s "$R/focal/dmel.fa" ]] && FOCAL_GENOME="$R/focal/dmel.fa"; } || true
    [[ -n "$FOCAL_GFF" ]]    || { [[ -s "$R/annotations/dmel.genes.gff" ]] && FOCAL_GFF="$R/annotations/dmel.genes.gff"; } || true
    [[ -n "$TE_DB" ]]        || { [[ -s "$R/te.fa" ]] && TE_DB="$R/te.fa"; } || true
    [[ -n "$NCRNA_DB" ]]     || { [[ -s "$R/ncrna.fa" ]] && NCRNA_DB="$R/ncrna.fa"; } || true
    log "Референсы: белки=${REF_PROT:-—} геномы=${REF_GEN:-—} HGT=${HGT_DB:-—}"
    log "           фокальный геном=${FOCAL_GENOME:-—} GFF=${FOCAL_GFF:-—} TE=${TE_DB:-—} ncRNA=${NCRNA_DB:-—}"
fi
[[ -z "$PROTEIN$LOCUS" || ( ${#GENES[@]} -le 1 && $ALL -eq 0 ) ]] \
    || die "--protein и --locus задаются только для одного гена"

activate_venv
require_cmd blastn blastp getorf makeblastdb
for f in extract_gene.py pick_orf.py classify_origin.py synteny.py $([[ $ONLINE -eq 1 ]] && echo ncbi_online.py); do
    [[ -f "$SCRIPTS_DIR/$f" ]] || die "Нет файла $SCRIPTS_DIR/$f — положите его рядом с classify_gene.sh"
done
if [[ $ONLINE -eq 1 ]]; then
    require_cmd tblastn curl
    GEN_SPECIES="${GEN_SPECIES:-${OUTGROUP:-dyak,dere,dana,dpse,dvir}}"
    GEN_TAXA="$(python3 "$SCRIPTS_DIR/ncbi_online.py" taxa "$GEN_SPECIES")" || exit 1
    HGT_TAXA='("Bacteria"[Organism] OR "Archaea"[Organism] OR "Viruses"[Organism] OR "Fungi"[Organism] OR "Viridiplantae"[Organism])'
    # маленькая библиотека мобильных элементов (Bergman lab, GitHub), один раз
    if [[ -z "$TE_DB" ]]; then
        TE_DB="$PROJECT_DIR/data/reference/te.fa"
        if [[ ! -s "$TE_DB" ]]; then
            log "Скачиваю библиотеку мобильных элементов (~1 МБ, один раз)..."
            mkdir -p "$PROJECT_DIR/data/reference/.te"
            if curl -fsSL --retry 3 "https://codeload.github.com/bergmanlab/drosophila-transposons/tar.gz/refs/heads/master" \
                    | tar -xz -C "$PROJECT_DIR/data/reference/.te"; then
                f="$(find "$PROJECT_DIR/data/reference/.te" -path '*current*' -iname 'transposon_sequence_set.fa*' | head -n1)"
                [[ -n "$f" ]] || f="$(find "$PROJECT_DIR/data/reference/.te" -iname '*transposon*.fa*' | sort | tail -n1)"
                [[ -n "$f" ]] && { case "$f" in *.gz) gzip -dc "$f" ;; *) cat "$f" ;; esac > "$TE_DB"; }
            fi
            rm -rf "$PROJECT_DIR/data/reference/.te"
            [[ -s "$TE_DB" ]] || { log "      не удалось — ветка TE будет недоступна"; TE_DB=""; }
        fi
    fi
    log "Онлайн-режим: белки в $ONLINE_SCOPE; геномы: $GEN_SPECIES"
fi
[[ -z "$REF_GEN" ]] || require_cmd tblastn
[[ -z "$FOCAL_GENOME" ]] || require_cmd blastdbcmd

CAND="$PROJECT_DIR/data/processed/$DATASET/candidates.fasta"
BLASTN_DB="$PROJECT_DIR/databases/blastn/$DATASET/cds"
BLASTP_DB="$PROJECT_DIR/databases/blastp/$DATASET/proteins"
REF_DB="$PROJECT_DIR/databases/reference"
[[ -s "$CAND" ]] || die "Датасет '$DATASET' не подготовлен. Запустите ./prepare_dataset.sh <файл>"

PROT_FMT="6 qseqid sseqid pident length qcovs evalue bitscore qlen slen qstart qend sstart send sseq"
NT_FMT="6 qseqid sseqid pident length qcovs evalue bitscore qlen slen qstart qend sstart send"
FLANK_FMT="6 qseqid sseqid pident length evalue bitscore sstart send"

regex_group() { python3 -c 'import re,sys; m=re.search(sys.argv[1],sys.argv[2]); print(m.group(1) if m else "")' "$1" "$2"; }

# ---------- Состав датасета по видам ----------
log "Виды в датасете '$DATASET':"
grep '^>' "$CAND" | sed 's/^>//; s/[[:space:]].*//' | python3 -c '
import re, sys, collections
rx = re.compile(sys.argv[1]); c = collections.Counter()
for l in sys.stdin:
    m = rx.search(l.strip()); c[m.group(1) if m else "?"] += 1
for k, v in c.most_common(): print(f"      {k}: {v}")' "$SPECIES_RE" >&2

# ---------- Базы (строятся один раз) ----------
fasta_in() { find "$1" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.faa' -o -name '*.fna' -o -name '*.fasta' \) | sort; }
stem() { local b; b="$(basename "$1")"; echo "${b%.*}"; }
build_db() {   # файл тип путь [доп. ключи]
    local f="$1" t="$2" db="$3"; shift 3
    local ext; [[ "$t" == prot ]] && ext=pin || ext=nin
    if [[ ! -s "$db.$ext" && ! -s "$db.00.$ext" ]]; then
        mkdir -p "$(dirname "$db")"
        log "Строю базу ($t): $(basename "$f")"
        makeblastdb -in "$f" -dbtype "$t" -out "$db" "$@" >/dev/null
    fi
}
for pair in "REF_PROT:prot:proteins" "REF_GEN:nucl:genomes" "HGT_DB:prot:hgt"; do
    IFS=: read -r var t sub <<< "$pair"
    dir="${!var}"; [[ -n "$dir" ]] || continue
    [[ -d "$dir" ]] || die "Нет каталога $dir"
    while read -r f; do build_db "$f" "$t" "$REF_DB/$sub/$(stem "$f")"; done < <(fasta_in "$dir")
done
[[ -z "$TE_DB" ]]    || { [[ -f "$TE_DB" ]] || die "Нет файла $TE_DB"; build_db "$TE_DB" nucl "$REF_DB/te/te"; }
[[ -z "$NCRNA_DB" ]] || { [[ -f "$NCRNA_DB" ]] || die "Нет файла $NCRNA_DB"; build_db "$NCRNA_DB" nucl "$REF_DB/ncrna/ncrna"; }
if [[ -n "$FOCAL_GENOME" ]]; then
    [[ -f "$FOCAL_GENOME" ]] || die "Нет файла $FOCAL_GENOME"
    build_db "$FOCAL_GENOME" nucl "$REF_DB/focal_genome/$(stem "$FOCAL_GENOME")" -parse_seqids
fi

# ---------- Один ген ----------
classify_one() {
    local GENE="$1" GS O P REC_ID FOC TASK TREE ORF_AA GENE_LEN
    GS="$(safe_name "$GENE")"
    O="$PROJECT_DIR/results/$DATASET/$GS/origin"; P="$O/$GS"
    mkdir -p "$O"; rm -f "$O"/*.tsv
    log "=== $GENE ==="

    REC_ID="$(python3 "$SCRIPTS_DIR/extract_gene.py" --fasta "$CAND" --query "$GENE" --out "$P.fasta")" \
        || { log "ОШИБКА: ген $GENE не найден в датасете"; return 1; }
    FOC="${FOCAL:-$(regex_group "$SPECIES_RE" "$REC_ID")}"
    [[ -n "$FOC" ]] || die "Не удалось определить вид из '$REC_ID' — задайте --focal"
    GENE_LEN="$(grep -v '^>' "$P.fasta" | tr -d '\n' | wc -c)"

    # [1] ORF
    log "[1/7] ORF (getorf)..."
    getorf -sequence "$P.fasta" -outseq "${P}_orfs_aa.fasta" -minsize 30 -find 1 2>/dev/null
    getorf -sequence "$P.fasta" -outseq "${P}_orfs_nt.fasta" -minsize 30 -find 3 2>/dev/null
    python3 "$SCRIPTS_DIR/pick_orf.py" "$P.fasta" "${P}_orfs_aa.fasta" "${P}_orfs_nt.fasta" "${P}_orf" \
        ${PROTEIN:+--protein "$PROTEIN"} && [[ -s "${P}_orf.tsv" ]] \
        || { log "ОШИБКА: не удалось определить ORF (getorf / pick_orf.py)"; return 1; }
    ORF_AA="$(awk -F'\t' 'NR==2{print $2}' "${P}_orf.tsv")"

    if [[ "$ORF_AA" -ge 1 ]]; then
        TASK=blastp; [[ "$ORF_AA" -lt 30 ]] && TASK=blastp-short
        # [2] Гомологи в датасете: белок (blastp против 6 рамок) и ДНК (blastn)
        log "[2/7] blastp / blastn против датасета..."
        blastp -task "$TASK" -query "${P}_orf.faa" -db "$BLASTP_DB" -evalue 10 -max_hsps 1 \
            -max_target_seqs "$MAX_TARGETS" -outfmt "$PROT_FMT" -num_threads "$THREADS" -out "${P}_cand_prot.tsv"
        blastn -task blastn -query "${P}_orf.fna" -db "$BLASTN_DB" -evalue 10 -max_hsps 1 \
            -max_target_seqs "$MAX_TARGETS" -outfmt "$NT_FMT" -num_threads "$THREADS" -out "${P}_cand_nt.tsv"

        # [3] Референсные протеомы и геномы
        log "[3/7] Референсы..."
        if [[ -n "$REF_PROT" ]]; then while read -r f; do local sp; sp="$(stem "$f")"
            blastp -task "$TASK" -query "${P}_orf.faa" -db "$REF_DB/proteins/$sp" -evalue 10 -max_hsps 1 \
                -max_target_seqs "$MAX_TARGETS" -outfmt "$PROT_FMT" -num_threads "$THREADS" \
                -out "${P}_refprot_${sp}.tsv"
        done < <(fasta_in "$REF_PROT"); fi
        if [[ -n "$REF_GEN" ]]; then while read -r f; do local sp; sp="$(stem "$f")"
            tblastn -query "${P}_orf.faa" -db "$REF_DB/genomes/$sp" -evalue 10 -max_hsps 1 \
                -max_target_seqs "$MAX_TARGETS" -outfmt "$PROT_FMT" -num_threads "$THREADS" \
                -out "${P}_prot_refgen_${sp}.tsv"
            blastn -task blastn -query "${P}_orf.fna" -db "$REF_DB/genomes/$sp" -evalue 10 -max_hsps 1 \
                -max_target_seqs "$MAX_TARGETS" -outfmt "$NT_FMT" -num_threads "$THREADS" \
                -out "${P}_nt_refgen_${sp}.tsv"
        done < <(fasta_in "$REF_GEN"); fi
        [[ -n "$REF_PROT$REF_GEN" ]] || log "      не заданы — только датасет"

        # [4] HGT: белки других ветвей жизни
        log "[4/7] Другие ветви жизни (HGT)..."
        if [[ -n "$HGT_DB" ]]; then while read -r f; do local g; g="$(stem "$f")"
            blastp -task "$TASK" -query "${P}_orf.faa" -db "$REF_DB/hgt/$g" -evalue 10 -max_hsps 1 \
                -max_target_seqs "$MAX_TARGETS" -outfmt "$PROT_FMT" -num_threads "$THREADS" \
                -out "${P}_hgt_${g}.tsv"
        done < <(fasta_in "$HGT_DB")
        elif [[ $ONLINE -eq 1 ]]; then log "      онлайн — см. [O3]"
        else log "      --hgt-db не задан — ветка 2.3.2 недоступна"; fi
    else
        log "[2-4/7] ORF не найден — поиск гомологов пропущен"
    fi

    # [O] Онлайн-поиск в NCBI (удалённый BLAST + API), с кэшем в origin/online/
    if [[ $ONLINE -eq 1 && "$ORF_AA" -ge 1 ]]; then
        local OD="$O/online" key gid loc NAME
        mkdir -p "$OD"
        remote() {   # remote <файл> <программа и аргументы...>: кэшируемый удалённый BLAST
            local out="$1"; shift
            if [[ -f "$out.done" && $REFRESH -eq 0 ]]; then log "      (из кэша) $(basename "$out")"; return 0; fi
            log "      удалённый $1 → $(basename "$out") (может занять несколько минут)..."
            if "$@" -remote -out "$out"; then touch "$out.done"
            else log "      ВНИМАНИЕ: удалённый поиск не удался ($(basename "$out"))"; : > "$out"; fi
        }
        NAME="$(python3 "$SCRIPTS_DIR/ncbi_online.py" name "$FOC")"
        log "[O1] NCBI Gene / Orthologs..."
        key="$(regex_group "$KEY_RE" "$REC_ID")"; gid=""
        if [[ -n "$key" ]]; then
            if [[ -s "$OD/gene_id.txt" && $REFRESH -eq 0 ]]; then gid="$(cat "$OD/gene_id.txt")"
            else gid="$(python3 "$SCRIPTS_DIR/ncbi_online.py" geneid "$key" "$NAME" || true)"
                 [[ -n "$gid" ]] && echo "$gid" > "$OD/gene_id.txt"; fi
        fi
        if [[ -n "$gid" ]]; then
            log "      GeneID: $gid"
            [[ -s "$OD/orthologs.tsv" && $REFRESH -eq 0 ]] \
                || python3 "$SCRIPTS_DIR/ncbi_online.py" orthologs "$gid" "$OD/orthologs.tsv" \
                || log "      ортологи не получены"
        else
            log "      ген не найден в NCBI Gene — без ортологов и синтении по фланкам"
        fi

        log "[O2] Белковые гомологи (refseq_protein, $ONLINE_SCOPE)..."
        remote "$OD/prot_scope.tsv" blastp -task "$TASK" -query "${P}_orf.faa" -db refseq_protein \
            -entrez_query "$ONLINE_SCOPE" -evalue 10 -max_target_seqs 200 -max_hsps 1 -outfmt "$PROT_FMT stitle"
        python3 "$SCRIPTS_DIR/ncbi_online.py" split "$OD/prot_scope.tsv" "$P" refprot

        log "[O3] Другие ветви жизни (HGT)..."
        remote "$OD/prot_hgt.tsv" blastp -task "$TASK" -query "${P}_orf.faa" -db refseq_protein \
            -entrez_query "$HGT_TAXA" -evalue 10 -max_target_seqs 50 -max_hsps 1 -outfmt "$PROT_FMT stitle"
        python3 "$SCRIPTS_DIR/ncbi_online.py" split "$OD/prot_hgt.tsv" "$P" hgt

        log "[O4] Геномы $GEN_SPECIES: tblastn / blastn (refseq_genomic)..."
        remote "$OD/genome_prot.tsv" tblastn -query "${P}_orf.faa" -db refseq_genomic \
            -entrez_query "$GEN_TAXA" -evalue 10 -max_target_seqs 50 -max_hsps 1 -outfmt "$PROT_FMT stitle"
        python3 "$SCRIPTS_DIR/ncbi_online.py" split "$OD/genome_prot.tsv" "$P" prot_refgen
        remote "$OD/genome_nt.tsv" blastn -task dc-megablast -query "${P}_orf.fna" -db refseq_genomic \
            -entrez_query "$GEN_TAXA" -evalue 10 -max_target_seqs 50 -max_hsps 1 -outfmt "$NT_FMT stitle"
        python3 "$SCRIPTS_DIR/ncbi_online.py" split "$OD/genome_nt.tsv" "$P" nt_refgen

        if [[ -n "$gid" ]]; then
            log "[O5] Синтения: фланки гена ($FLANK нт) в геномах $GEN_SPECIES..."
            if [[ ! -s "$OD/flanks.fna" || $REFRESH -eq 1 ]]; then
                loc="$(python3 "$SCRIPTS_DIR/ncbi_online.py" locus "$gid" || true)"
                if [[ -n "$loc" ]]; then
                    local acc s e; IFS=$'\t' read -r acc s e <<< "$loc"
                    log "      локус: $acc:$s-$e"; echo "$loc" > "$OD/locus.tsv"
                    python3 "$SCRIPTS_DIR/ncbi_online.py" flanks "$acc" "$s" "$e" "$FLANK" "$OD/flanks.fna" || true
                fi
            fi
            if [[ -s "$OD/flanks.fna" ]]; then
                remote "$OD/flanks_hits.tsv" blastn -task dc-megablast -query "$OD/flanks.fna" -db refseq_genomic \
                    -entrez_query "$GEN_TAXA" -evalue 1e-10 -max_target_seqs 20 -max_hsps 5 -outfmt "$FLANK_FMT stitle"
                python3 "$SCRIPTS_DIR/ncbi_online.py" split "$OD/flanks_hits.tsv" "$P" flanks
                local glen; IFS=$'\t' read -r _ s e < "$OD/locus.tsv"; glen=$(( e - s + 1 ))
                for ff in "$O"/*_flanks_*.tsv; do
                    [[ -e "$ff" ]] || continue
                    local sp="${ff##*_flanks_}"; sp="${sp%.tsv}"
                    python3 "$SCRIPTS_DIR/synteny.py" window "$ff" "$sp" "$glen" "$FLANK" "$O/synteny_${sp}.tsv"
                done
            fi
        fi

        log "[O6] ncRNA $NAME (refseq_rna)..."
        remote "$OD/ncrna.tsv" blastn -task blastn -query "$P.fasta" -db refseq_rna \
            -entrez_query "\"$NAME\"[Organism] AND biomol_ncrna[PROP]" -evalue 1e-5 \
            -max_target_seqs 20 -max_hsps 1 -outfmt "$NT_FMT stitle"
        python3 "$SCRIPTS_DIR/ncbi_online.py" split "$OD/ncrna.tsv" "$P" ncrna
    fi

    # [5] Синтения через фланки: геном фокального вида + координаты гена → окна у видов из --ref-genomes
    log "[5/7] Синтения..."
    if [[ $ONLINE -eq 1 ]]; then
        log "      онлайн: окна по фланкам и ортологи NCBI (см. [O1], [O5])"
    elif [[ -n "$FOCAL_GENOME" && -n "$REF_GEN" && ( -n "$LOCUS" || -n "$FOCAL_GFF" ) ]]; then
        local chr s e key fdb loc
        if [[ -n "$LOCUS" ]]; then
            [[ "$LOCUS" =~ ^(.+):([0-9]+)-([0-9]+)$ ]] || die "--locus: формат CHR:START-END"
            chr="${BASH_REMATCH[1]}"; s="${BASH_REMATCH[2]}"; e="${BASH_REMATCH[3]}"
        else
            key="$(regex_group "$KEY_RE" "$REC_ID")"
            [[ -n "$key" ]] || die "В '$REC_ID' нет ключа по --gene-key-regex — задайте --locus"
            loc="$(python3 "$SCRIPTS_DIR/synteny.py" locus "$FOCAL_GFF" "$key")"
            IFS=$'\t' read -r chr s e _ <<< "$loc"
        fi
        log "      локус: $chr:$s-$e, фланки по $FLANK нт"
        fdb="$REF_DB/focal_genome/$(stem "$FOCAL_GENOME")"
        { blastdbcmd -db "$fdb" -entry "$chr" -range "$(( s > FLANK ? s - FLANK : 1 ))-$(( s > 1 ? s - 1 : 1 ))" \
              | sed '1s/.*/>left/'
          blastdbcmd -db "$fdb" -entry "$chr" -range "$(( e + 1 ))-$(( e + FLANK ))" \
              | sed '1s/.*/>right/'; } > "${P}_flanks.fna" \
            || log "      не удалось извлечь фланки (проверьте имя хромосомы)"
        while read -r f; do local sp; sp="$(stem "$f")"
            blastn -task dc-megablast -query "${P}_flanks.fna" -db "$REF_DB/genomes/$sp" -evalue 1e-10 \
                -max_target_seqs 5 -max_hsps 5 -outfmt "$FLANK_FMT" -num_threads "$THREADS" \
                -out "${P}_flanks_${sp}.tsv"
            python3 "$SCRIPTS_DIR/synteny.py" window "${P}_flanks_${sp}.tsv" "$sp" \
                "$(( e - s + 1 ))" "$FLANK" "$O/synteny_${sp}.tsv"
        done < <(fasta_in "$REF_GEN")
    else
        log "      по фланкам — нет (нужны --focal-genome, --ref-genomes и --focal-gff/--locus);"
        log "      в датасете синтенными считаются записи с тем же ключом гена ($KEY_RE)"
    fi

    # [6] TE / ncRNA: blastn всего гена
    log "[6/7] Связь с TE / ncRNA..."
    [[ -z "$TE_DB" ]] || blastn -task blastn -query "$P.fasta" -db "$REF_DB/te/te" -evalue 1e-5 \
        -max_target_seqs 50 -max_hsps 1 -outfmt "$NT_FMT" -num_threads "$THREADS" -out "${P}_te.tsv"
    [[ -z "$NCRNA_DB" ]] || blastn -task blastn -query "$P.fasta" -db "$REF_DB/ncrna/ncrna" -evalue 1e-5 \
        -max_target_seqs 50 -max_hsps 1 -outfmt "$NT_FMT" -num_threads "$THREADS" -out "${P}_ncrna.tsv"
    [[ -n "$TE_DB$NCRNA_DB" || $ONLINE -eq 1 ]] || log "      --te-db / --ncrna-db не заданы — ветка 3.2 недоступна"

    # [7] Дерево (для спорных случаев дупликация/древний)
    TREE="$PROJECT_DIR/results/$DATASET/$GS/phylogeny/${GS}_blastp.treefile"
    if [[ ! -s "$TREE" && $BUILD_TREE -eq 1 ]]; then
        log "[7/7] Строю дерево (run_gene.sh)..."
        "$PROJECT_DIR/run_gene.sh" -d "$DATASET" -g "$GENE" -t "$THREADS" || log "      не удалось"
    fi
    [[ -s "$TREE" ]] && log "[7/7] Дерево: есть" || { log "[7/7] Дерева нет (--tree)"; TREE=""; }

    # Классификация
    shopt -s nullglob
    local prot=("$O"/*_cand_prot.tsv "$O"/*_refprot_*.tsv "$O"/*_prot_refgen_*.tsv)
    local nt=("$O"/*_cand_nt.tsv "$O"/*_nt_refgen_*.tsv)
    local hgt=("$O"/*_hgt_*.tsv) te=("$O"/*_te.tsv) nc=("$O"/*_ncrna.tsv) syn=("$O"/synteny_*.tsv)
    shopt -u nullglob
    local RESULT
    RESULT="$(python3 "$SCRIPTS_DIR/classify_origin.py" \
        --rec-id "$REC_ID" --focal "$FOC" --ingroup "$INGROUP" \
        --outgroup "$([[ $ONLINE -eq 1 ]] || echo "$OUTGROUP")" \
        --orthologs "$O/online/orthologs.tsv" \
        --species-regex "$SPECIES_RE" --gene-key-regex "$KEY_RE" --orf-tsv "${P}_orf.tsv" \
        --prot-hits "${prot[@]}" --nt-hits "${nt[@]}" --hgt-hits "${hgt[@]}" \
        --te-hits "${te[@]}" --ncrna-hits "${nc[@]}" --synteny "${syn[@]}" \
        ${TREE:+--treefile "$TREE"} --outdir "$O" $([[ $REQ_SYN -eq 1 ]] && echo --require-synteny) \
        --evalue "$EVALUE" --min-pident "$MIN_PIDENT" --min-cov "$MIN_COV" --part-cov "$PART_COV" \
        --max-overlap "$MAX_OVERLAP" --min-orf-aa "$MIN_ORF_AA" --nt-pident "$NT_PIDENT" \
        --nt-cov "$NT_COV" --min-ai "$MIN_AI" --te-evalue "$TE_EVALUE" --te-pident "$TE_PIDENT" \
        --te-cov "$TE_COV" --min-uf "$MIN_UF" --min-alrt "$MIN_ALRT")" \
        || { log "ОШИБКА: классификатор завершился с ошибкой (см. выше)"; return 1; }

    local SUM="$PROJECT_DIR/results/$DATASET/origin_summary.tsv"
    [[ -s "$SUM" ]] || printf 'gene\tverdict\tconfidence\treason\n' > "$SUM"
    awk -F'\t' -v g="$REC_ID" '$1 != g' "$SUM" > "$SUM.tmp" && mv "$SUM.tmp" "$SUM"
    echo "$RESULT" >> "$SUM"
    log "ИТОГ: $(awk -F'\t' '{print $2 " (уверенность: " $3 ")"}' <<< "$RESULT")"
    log "Причина: $(cut -f4 <<< "$RESULT")"
    log "Отчёт: $O/origin_report.md"
}

[[ $ALL -eq 0 ]] || mapfile -t GENES < <(grep '^>' "$CAND" | sed 's/^>//; s/[[:space:]].*//')

n_fail=0
for g in "${GENES[@]}"; do
    ( classify_one "$g" ) || { n_fail=$((n_fail + 1)); log "Ошибка на гене $g"; }
done
[[ ${#GENES[@]} -le 1 ]] || log "Сводка: $PROJECT_DIR/results/$DATASET/origin_summary.tsv"
[[ $n_fail -eq 0 ]]
