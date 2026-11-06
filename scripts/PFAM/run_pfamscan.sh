#!/usr/bin/env bash

set -euo pipefail

# Usage:
#   run_pfamscan.sh INPUT_DIR OUTPUT_DIR PFAM_DB THREADS
#       Processes all FASTA files in INPUT_DIR.
#
#   run_pfamscan.sh INPUT_DIR OUTPUT_DIR PFAM_DB THREADS file1.fa file2.fa ...
#       Processes only the explicitly supplied FASTA files.
#
# The OrthoSimulator PFAM wrapper normally calls pfam_scan.pl directly, because it must
# process only the current --Proteome input. This helper is kept for manual use
# and mirrors the same current-input-only behaviour when FASTA files are given.

INPUT_DIR="${1:-pfam_genomes}"
OUTPUT_DIR="${2:-pfamscan_results}"
PFAM_DB="${3:-pfam_db}"
THREADS="${4:-1}"
shift 4 || true

mkdir -p "$OUTPUT_DIR"

FASTAS=()

if [ "$#" -gt 0 ]; then
    for fasta in "$@"; do
        FASTAS+=("$fasta")
    done
else
    shopt -s nullglob
    FASTAS=(
        "$INPUT_DIR"/*.fa
        "$INPUT_DIR"/*.fasta
        "$INPUT_DIR"/*.faa
        "$INPUT_DIR"/*.fas
        "$INPUT_DIR"/*.fna
    )
fi

if [ ${#FASTAS[@]} -eq 0 ]; then
    echo "ERROR: no FASTA files found or supplied" >&2
    exit 1
fi

for fasta in "${FASTAS[@]}"; do
    if [ ! -f "$fasta" ]; then
        echo "ERROR: FASTA file not found: $fasta" >&2
        exit 1
    fi

    filename=$(basename "$fasta")
    stem="${filename%.*}"
    outfile="$OUTPUT_DIR/${stem}_pfam.txt"

    echo "Processing $fasta -> $outfile"

    pfam_scan.pl \
        -fasta "$fasta" \
        -dir "$PFAM_DB" \
        -cpu "$THREADS" \
        -outfile "$outfile"
done

echo "All requested files processed."
