#!/usr/bin/env bash

set -euo pipefail

PFAM_URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
OUTDIR="${1:-pfam_db}"

mkdir -p "$OUTDIR"

echo "Downloading Pfam files to: $OUTDIR"
wget -O "$OUTDIR/Pfam-A.hmm.dat.gz" "$PFAM_URL/Pfam-A.hmm.dat.gz"
wget -O "$OUTDIR/Pfam-A.hmm.gz" "$PFAM_URL/Pfam-A.hmm.gz"

echo "Unpacking files..."
gunzip -c "$OUTDIR/Pfam-A.hmm.dat.gz" > "$OUTDIR/Pfam-A.hmm.dat"
gunzip -c "$OUTDIR/Pfam-A.hmm.gz" > "$OUTDIR/Pfam-A.hmm"

echo "Removing compressed files..."
rm -f "$OUTDIR/Pfam-A.hmm.dat.gz" "$OUTDIR/Pfam-A.hmm.gz"

echo "Preparing Pfam database for HMMER..."
hmmpress -f "$OUTDIR/Pfam-A.hmm"

echo "Done. Pfam database is ready in: $OUTDIR"
