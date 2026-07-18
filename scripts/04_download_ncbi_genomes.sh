#!/usr/bin/env bash

set -euo pipefail

ACCESSION_FILE=""
OUTPUT_DIR=""
INCLUDE="genome,protein,gff3"
DEHYDRATED=false
FORCE=false

usage() {
    cat <<'EOF'
Usage:
  bash scripts/04_download_ncbi_genomes.sh \
      --accessions ACCESSIONS.txt \
      --output-dir OUTPUT_DIRECTORY

Options:
  --accessions FILE
      Text file containing one GCA_/GCF_ accession per line.

  --output-dir DIRECTORY
      Output directory for the NCBI Datasets package.

  --include TYPES
      Comma-separated data types.
      Default: genome,protein,gff3

  --dehydrated
      Create a dehydrated package and then rehydrate it.

  --force
      Remove an existing output directory before downloading.

  -h, --help
      Show this help message.
EOF
}


while [[ $# -gt 0 ]]; do
    case "$1" in
        --accessions)
            ACCESSION_FILE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --include)
            INCLUDE="$2"
            shift 2
            ;;
        --dehydrated)
            DEHYDRATED=true
            shift
            ;;
        --force)
            FORCE=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done


if [[ -z "$ACCESSION_FILE" ]]; then
    echo "[ERROR] --accessions is required." >&2
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "[ERROR] --output-dir is required." >&2
    exit 1
fi

if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "[ERROR] Accession file not found: $ACCESSION_FILE" >&2
    exit 1
fi

if ! command -v datasets >/dev/null 2>&1; then
    echo "[ERROR] NCBI datasets command was not found." >&2
    echo "Install with:" >&2
    echo "  mamba install -c conda-forge ncbi-datasets-cli" >&2
    exit 1
fi

if ! command -v unzip >/dev/null 2>&1; then
    echo "[ERROR] unzip command was not found." >&2
    exit 1
fi


ACCESSION_COUNT=$(
    grep -Ev '^[[:space:]]*$' "$ACCESSION_FILE" |
    sort -u |
    wc -l |
    tr -d ' '
)

if [[ "$ACCESSION_COUNT" -eq 0 ]]; then
    echo "[ERROR] Accession file is empty." >&2
    exit 1
fi


if [[ -e "$OUTPUT_DIR" ]]; then
    if [[ "$FORCE" == true ]]; then
        echo "[REMOVE] Existing output: $OUTPUT_DIR"
        rm -rf "$OUTPUT_DIR"
    else
        echo "[ERROR] Output already exists: $OUTPUT_DIR" >&2
        echo "Use --force to replace it." >&2
        exit 1
    fi
fi


mkdir -p "$OUTPUT_DIR"

ZIP_PATH="${OUTPUT_DIR}/ncbi_dataset.zip"
EXTRACT_DIR="${OUTPUT_DIR}/package"
LOG_PATH="${OUTPUT_DIR}/download.log"


echo "========================================================================"
echo "GutSporePredict v4.0-alpha1"
echo "NCBI genome download"
echo "========================================================================"
echo "Accession file:  $ACCESSION_FILE"
echo "Accessions:      $ACCESSION_COUNT"
echo "Output:          $OUTPUT_DIR"
echo "Included files:  $INCLUDE"
echo "Dehydrated:      $DEHYDRATED"
echo


DOWNLOAD_COMMAND=(
    datasets
    download
    genome
    accession
    --inputfile "$ACCESSION_FILE"
    --include "$INCLUDE"
    --filename "$ZIP_PATH"
)


if [[ "$DEHYDRATED" == true ]]; then
    DOWNLOAD_COMMAND+=(--dehydrated)
fi


{
    echo "[COMMAND] ${DOWNLOAD_COMMAND[*]}"
    "${DOWNLOAD_COMMAND[@]}"
} 2>&1 | tee "$LOG_PATH"


echo
echo "[ZIP TEST] $ZIP_PATH"
unzip -t "$ZIP_PATH" >/dev/null

mkdir -p "$EXTRACT_DIR"

echo "[EXTRACT] $EXTRACT_DIR"
unzip -q "$ZIP_PATH" -d "$EXTRACT_DIR"


if [[ "$DEHYDRATED" == true ]]; then
    echo "[REHYDRATE] $EXTRACT_DIR"

    datasets rehydrate \
        --directory "$EXTRACT_DIR" \
        2>&1 | tee -a "$LOG_PATH"
fi


DATA_DIR="${EXTRACT_DIR}/ncbi_dataset/data"

if [[ ! -d "$DATA_DIR" ]]; then
    echo "[ERROR] NCBI data directory was not created:" >&2
    echo "  $DATA_DIR" >&2
    exit 1
fi


GENOME_COUNT=$(
    find "$DATA_DIR" \
        -type f \
        -name '*_genomic.fna' |
    wc -l |
    tr -d ' '
)

PROTEIN_COUNT=$(
    find "$DATA_DIR" \
        -type f \
        -name 'protein.faa' |
    wc -l |
    tr -d ' '
)

GFF_COUNT=$(
    find "$DATA_DIR" \
        -type f \
        -name 'genomic.gff' |
    wc -l |
    tr -d ' '
)


echo
echo "========================================================================"
echo "[SUCCESS] Download package completed."
echo "Requested assemblies: $ACCESSION_COUNT"
echo "Genome FASTA files:   $GENOME_COUNT"
echo "Protein FASTA files:  $PROTEIN_COUNT"
echo "GFF3 files:           $GFF_COUNT"
echo "Package directory:    $EXTRACT_DIR"
echo "Log:                  $LOG_PATH"
echo "========================================================================"
