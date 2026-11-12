#!/usr/bin/env bash
# =============================================================================
# fetch_references.sh — скачивает референсные данные для classify_gene.sh
#
#   FlyBase          геном, гены (GFF), белки, ncRNA, мобильные элементы D. melanogaster
#   NCBI Datasets    геномы и белки других видов Drosophila (RefSeq reference)
#   UniProt          Swiss-Prot белки других ветвей жизни (для проверки HGT)
#
# Использование:
#   ./fetch_references.sh [опции]
# Опции:
#   --species "dsim dyak ..."   виды для NCBI (по умолчанию: dsim dsec dyak dere dana dpse dvir)
#   --no-flybase                не качать данные D. melanogaster (FlyBase, запасной — NCBI)
#   --no-ncbi                   не качать геномы других видов
#   --no-hgt                    не качать UniProt
#   --hgt-groups "bacteria ..." группы для HGT (по умолчанию: bacteria archaea viruses fungi plants)
#   --force                     перекачать уже скачанное
#
# Результат (data/reference/):
#   focal/dmel.fa                 геном D. melanogaster        → --focal-genome
#   annotations/dmel.genes.gff    только строки gene из GFF   → --focal-gff
#   proteins/<вид>.faa            протеомы (dmel + NCBI)       → --ref-proteins
#   genomes/<вид>.fa              геномы других видов          → --ref-genomes
#   hgt/<группа>.faa              Swiss-Prot по группам        → --hgt-db
#   te.fa, ncrna.fa               мобильные элементы, ncRNA    → --te-db, --ncrna-db
#   SOURCES.tsv                   откуда и когда скачано (для отчёта)
# После загрузки: ./classify_gene.sh ... --use-references
#
# Объём: геномы ~150–250 МБ на вид, FlyBase ~300 МБ, UniProt ~200 МБ.
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

usage() { sed -n '2,33p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

SPECIES="dsim dsec dyak dere dana dpse dvir"
HGT_GROUPS="bacteria archaea viruses fungi plants"
DO_FB=1; DO_NCBI=1; DO_HGT=1; FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --species)     SPECIES="${2:?}"; shift 2 ;;
        --hgt-groups)  HGT_GROUPS="${2:?}"; shift 2 ;;
        --no-flybase)  DO_FB=0; shift ;;
        --no-ncbi)     DO_NCBI=0; shift ;;
        --no-hgt)      DO_HGT=0; shift ;;
        --force)       FORCE=1; shift ;;
        -h|--help)     usage; exit 0 ;;
        *)             die "Неизвестный аргумент: $1" ;;
    esac
done

require_cmd curl gzip unzip
REF="$PROJECT_DIR/data/reference"
TMP="$PROJECT_DIR/data/reference/.tmp"
mkdir -p "$REF"/{focal,annotations,proteins,genomes,hgt} "$TMP"
SOURCES="$REF/SOURCES.tsv"
[[ -s "$SOURCES" ]] || printf 'file\tsource\turl_or_query\tdate\n' > "$SOURCES"

record() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$(date +%F)" >> "$SOURCES"; }
need() { [[ $FORCE -eq 1 || ! -s "$1" ]]; }
fetch() { curl -fL --retry 3 --retry-delay 5 -C - --progress-bar -o "$2" "$1"; }

declare -A TAXON=(
    [dsim]="Drosophila simulans"     [dsec]="Drosophila sechellia"
    [dyak]="Drosophila yakuba"       [dere]="Drosophila erecta"
    [dana]="Drosophila ananassae"    [dpse]="Drosophila pseudoobscura"
    [dper]="Drosophila persimilis"   [dwil]="Drosophila willistoni"
    [dmoj]="Drosophila mojavensis"   [dvir]="Drosophila virilis"
    [dgri]="Drosophila grimshawi"
)
declare -A HGT_TAXID=( [bacteria]=2 [archaea]=2157 [viruses]=10239 [fungi]=4751 [plants]=33090 )

# ============================== FlyBase ==============================
# Пробуем несколько адресов FlyBase; если ни один не отвечает — берём dmel из NCBI.
install_datasets() {
    [[ -n "${DATASETS:-}" ]] && return 0
    DATASETS="$(command -v datasets || true)"
    [[ -z "$DATASETS" && -x "$VENV_DIR/bin/datasets" ]] && DATASETS="$VENV_DIR/bin/datasets"
    if [[ -z "$DATASETS" ]]; then
        log "Устанавливаю NCBI datasets CLI в .venv/bin..."
        mkdir -p "$VENV_DIR/bin"
        fetch "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets" "$VENV_DIR/bin/datasets"
        chmod +x "$VENV_DIR/bin/datasets"; DATASETS="$VENV_DIR/bin/datasets"
    fi
}

fetch_flybase() {
    local FB="" LIST="" base
    for base in https://ftp.flybase.net https://ftp.flybase.org https://s3ftp.flybase.org; do
        log "FlyBase: пробую $base ..."
        if LIST="$(curl -fsSL --connect-timeout 20 --max-time 120 \
                   "$base/genomes/Drosophila_melanogaster/current/fasta/" 2>&1)" \
           && grep -q 'dmel-all-chromosome' <<< "$LIST"; then
            FB="$base/genomes/Drosophila_melanogaster/current"; break
        fi
        log "      не отвечает: $(head -c 200 <<< "$LIST" | tr '\n' ' ')"
    done
    [[ -n "$FB" ]] || return 1

    fb_file() { grep -oE "$1" <<< "$LIST" | sort -u | head -n1; }
    local CHR TRN NCR TEF GFF
    CHR="$(fb_file 'dmel-all-chromosome-r[0-9.]+\.fasta\.gz')"
    TRN="$(fb_file 'dmel-all-translation-r[0-9.]+\.fasta\.gz')"
    NCR="$(fb_file 'dmel-all-ncRNA-r[0-9.]+\.fasta\.gz')"
    TEF="$(fb_file 'dmel-all-transposon-r[0-9.]+\.fasta\.gz')"
    GFF="$(curl -fsL --connect-timeout 20 "$FB/gff/" | grep -oE 'dmel-all-r[0-9.]+\.gff\.gz' | sort -u | head -n1 || true)"
    log "Релиз: ${CHR:-?} | ${GFF:-?}"

    if [[ -n "$CHR" ]] && need "$REF/focal/dmel.fa"; then
        log "Геном dmel..."
        fetch "$FB/fasta/$CHR" "$TMP/$CHR"
        gzip -dc "$TMP/$CHR" > "$REF/focal/dmel.fa"; record focal/dmel.fa FlyBase "$FB/fasta/$CHR"
    fi
    if [[ -n "$GFF" ]] && need "$REF/annotations/dmel.genes.gff"; then
        log "GFF dmel (оставляю только строки gene)..."
        curl -fL --retry 3 --progress-bar "$FB/gff/$GFF" | gzip -dc \
            | awk -F'\t' 'NF >= 9 && $3 == "gene"' > "$REF/annotations/dmel.genes.gff"
        record annotations/dmel.genes.gff FlyBase "$FB/gff/$GFF"
    fi
    if [[ -n "$TRN" ]] && need "$REF/proteins/dmel.faa"; then
        log "Белки dmel..."
        fetch "$FB/fasta/$TRN" "$TMP/$TRN"
        # >FBpp... parent=FBgn...,FBtr... → >FBgn..._FBpp... (FBgn нужен, чтобы узнать сам ген)
        gzip -dc "$TMP/$TRN" | awk '/^>/ { id = substr($1, 2); g = ""
                if (match($0, /parent=FBgn[0-9]+/)) g = substr($0, RSTART + 7, RLENGTH - 7)
                print ">" (g != "" ? g "_" : "") id; next } { print }' > "$REF/proteins/dmel.faa"
        record proteins/dmel.faa FlyBase "$FB/fasta/$TRN"
    fi
    if [[ -n "$NCR" ]] && need "$REF/ncrna.fa"; then
        log "ncRNA dmel..."
        fetch "$FB/fasta/$NCR" "$TMP/$NCR"; gzip -dc "$TMP/$NCR" > "$REF/ncrna.fa"
        record ncrna.fa FlyBase "$FB/fasta/$NCR"
    fi
    if [[ -n "$TEF" ]] && need "$REF/te.fa"; then
        log "Мобильные элементы dmel..."
        fetch "$FB/fasta/$TEF" "$TMP/$TEF"; gzip -dc "$TMP/$TEF" > "$REF/te.fa"
        record te.fa FlyBase "$FB/fasta/$TEF"
    fi
    return 0
}

# Запасной вариант: D. melanogaster из NCBI (геном, GFF, белки, ncRNA — без мобильных элементов)
fetch_dmel_ncbi() {
    if ! need "$REF/focal/dmel.fa" && [[ -s "$REF/proteins/dmel.faa" && -s "$REF/annotations/dmel.genes.gff" ]]; then
        log "D. melanogaster уже скачан"; return 0
    fi
    install_datasets
    log "NCBI: Drosophila melanogaster (геном, GFF, белки, РНК)..."
    rm -rf "$TMP/dmel" "$TMP/dmel.zip"
    "$DATASETS" download genome taxon "Drosophila melanogaster" --reference \
        --include genome,protein,gff3,rna --filename "$TMP/dmel.zip" || return 1
    unzip -q -o "$TMP/dmel.zip" -d "$TMP/dmel"
    local gen gff prot rna acc
    gen="$(find "$TMP/dmel" -name '*_genomic.fna' | head -n1)"
    gff="$(find "$TMP/dmel" -name 'genomic.gff' | head -n1)"
    prot="$(find "$TMP/dmel" -name 'protein.faa' | head -n1)"
    rna="$(find "$TMP/dmel" -name 'rna.fna' | head -n1)"
    acc="$(basename "$(dirname "${gen:-x/x}")")"
    [[ -n "$gen" ]]  && { mv "$gen" "$REF/focal/dmel.fa"; record focal/dmel.fa "NCBI $acc" "D. melanogaster"; }
    # в GFF NCBI ключ FlyBase лежит в атрибуте Dbxref=FLYBASE:FBgn... — synteny.py его находит
    [[ -n "$gff" ]]  && { awk -F'\t' 'NF >= 9 && $3 == "gene"' "$gff" > "$REF/annotations/dmel.genes.gff"
                          record annotations/dmel.genes.gff "NCBI $acc" "genomic.gff (gene)"; }
    [[ -n "$prot" ]] && { mv "$prot" "$REF/proteins/dmel.faa"; record proteins/dmel.faa "NCBI $acc" "protein.faa"; }
    # некодирующие РНК RefSeq имеют префиксы NR_ / XR_
    [[ -n "$rna" ]]  && { awk '/^>/ { keep = ($1 ~ /^>(NR|XR)_/) } keep' "$rna" > "$REF/ncrna.fa"
                          record ncrna.fa "NCBI $acc" "rna.fna (NR_/XR_)"; }
    rm -rf "$TMP/dmel" "$TMP/dmel.zip"
}

# Мобильные элементы без FlyBase: канонический набор транспозонов Drosophila (Bergman lab, GitHub)
fetch_te_github() {
    need "$REF/te.fa" || { log "te.fa уже есть"; return 0; }
    local url="https://codeload.github.com/bergmanlab/drosophila-transposons/tar.gz/refs/heads/master" f
    log "Мобильные элементы: bergmanlab/drosophila-transposons (GitHub)..."
    rm -rf "$TMP/te"; mkdir -p "$TMP/te"
    curl -fsSL --retry 3 "$url" | tar -xz -C "$TMP/te" || { log "ВНИМАНИЕ: не удалось скачать с GitHub"; return 1; }
    # весь род (transposon_sequence_set.fa), иначе только D. melanogaster
    f="$(find "$TMP/te" -path '*current*' -iname 'transposon_sequence_set.fa*' | head -n1)"
    [[ -n "$f" ]] || f="$(find "$TMP/te" -path '*current*' -iname '*transposon*.fa*' | head -n1)"
    [[ -n "$f" ]] || f="$(find "$TMP/te" -iname '*transposon*.fa*' | sort | tail -n1)"
    [[ -n "$f" ]] || { log "ВНИМАНИЕ: в репозитории не найден FASTA с транспозонами"; return 1; }
    case "$f" in *.gz) gzip -dc "$f" ;; *) cat "$f" ;; esac > "$REF/te.fa"
    record te.fa "bergmanlab/drosophila-transposons" "${f#$TMP/te/}"
    log "      элементов: $(grep -c '^>' "$REF/te.fa")"
    rm -rf "$TMP/te"
}

if [[ $DO_FB -eq 1 ]]; then
    if ! fetch_flybase; then
        log "FlyBase недоступен ни по одному адресу — беру D. melanogaster из NCBI."
        fetch_dmel_ncbi || log "ВНИМАНИЕ: не удалось скачать D. melanogaster и из NCBI"
    fi
    [[ -s "$REF/te.fa" ]] || fetch_te_github || true
fi

# ============================== NCBI ==============================
if [[ $DO_NCBI -eq 1 ]]; then
    install_datasets
    for sp in $SPECIES; do
        name="${TAXON[$sp]:-}"
        [[ -n "$name" ]] || { log "Неизвестный вид '$sp' (есть: ${!TAXON[*]}) — пропускаю"; continue; }
        if ! need "$REF/genomes/$sp.fa" && [[ -s "$REF/proteins/$sp.faa" ]]; then
            log "$sp уже скачан"; continue
        fi
        log "NCBI: $name ($sp)..."
        rm -rf "$TMP/$sp" "$TMP/$sp.zip"
        if ! "$DATASETS" download genome taxon "$name" --reference --include genome,protein \
                --filename "$TMP/$sp.zip"; then
            log "ВНИМАНИЕ: не удалось скачать $name"; continue
        fi
        unzip -q -o "$TMP/$sp.zip" -d "$TMP/$sp"
        gen="$(find "$TMP/$sp" -name '*_genomic.fna' | head -n1)"
        prot="$(find "$TMP/$sp" -name 'protein.faa' | head -n1)"
        acc="$(basename "$(dirname "${gen:-x/x}")")"
        if [[ -n "$gen" ]]; then mv "$gen" "$REF/genomes/$sp.fa"; record "genomes/$sp.fa" "NCBI $acc" "$name"; fi
        if [[ -n "$prot" ]]; then mv "$prot" "$REF/proteins/$sp.faa"; record "proteins/$sp.faa" "NCBI $acc" "$name"
        else log "      у сборки $name нет аннотации белков — только геном"; fi
        rm -rf "$TMP/$sp" "$TMP/$sp.zip"
    done
fi

# ============================== UniProt (HGT) ==============================
uniprot_rest() {   # $1 группа $2 taxid → FASTA в $REF/hgt/$1.faa
    local q="(taxonomy_id:$2) AND (reviewed:true)" i
    for i in 1 2 3; do
        log "      REST, попытка $i..."
        if curl -fL --http1.1 --retry 2 --progress-bar -G "https://rest.uniprot.org/uniprotkb/stream" \
                --data-urlencode "query=$q" --data-urlencode "format=fasta" \
                --data-urlencode "compressed=true" -o "$TMP/$1.faa.gz" \
           && gzip -t "$TMP/$1.faa.gz" 2>/dev/null; then
            gzip -dc "$TMP/$1.faa.gz" > "$REF/hgt/$1.faa"; rm -f "$TMP/$1.faa.gz"
            record "hgt/$1.faa" "UniProt Swiss-Prot REST" "$q"; return 0
        fi
        sleep 5
    done
    return 1
}

uniprot_ftp() {    # $1 группа = имя таксономического раздела Swiss-Prot
    local m path="databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_$1.dat.gz"
    for m in https://ftp.uniprot.org/pub https://ftp.ebi.ac.uk/pub https://ftp.expasy.org; do
        log "      FTP-зеркало $m ..."
        if fetch "$m/$path" "$TMP/$1.dat.gz" && gzip -t "$TMP/$1.dat.gz" 2>/dev/null; then
            # Swiss-Prot flat file → FASTA: ID-строка → заголовок, строки с отступом → последовательность
            gzip -dc "$TMP/$1.dat.gz" | awk '/^ID / { print ">" $2; next }
                                             /^     / { gsub(/ /, ""); print }' > "$REF/hgt/$1.faa"
            rm -f "$TMP/$1.dat.gz"
            record "hgt/$1.faa" "UniProt Swiss-Prot FTP" "$m/$path"; return 0
        fi
        rm -f "$TMP/$1.dat.gz"
    done
    return 1
}

if [[ $DO_HGT -eq 1 ]]; then
    for g in $HGT_GROUPS; do
        tax="${HGT_TAXID[$g]:-}"
        [[ -n "$tax" ]] || { log "Неизвестная группа '$g' (есть: ${!HGT_TAXID[*]})"; continue; }
        need "$REF/hgt/$g.faa" || { log "HGT $g уже скачан"; continue; }
        log "UniProt Swiss-Prot: $g (taxonomy_id:$tax)..."
        if uniprot_rest "$g" "$tax" || uniprot_ftp "$g"; then
            log "      последовательностей: $(grep -c '^>' "$REF/hgt/$g.faa")"
        else
            rm -f "$REF/hgt/$g.faa"
            log "ВНИМАНИЕ: $g не скачан ни через REST, ни через FTP — пропускаю"
        fi
    done
fi

rm -rf "$TMP"
log "Готово. Содержимое $REF:"
du -sh "$REF"/* 2>/dev/null | sed 's/^/      /' >&2
log "Запуск анализа с референсами:  ./classify_gene.sh -d ДАТАСЕТ -g ГЕН --use-references"
