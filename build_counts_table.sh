#!/bin/bash
set -euo pipefail

# ===== Define Directories =====
BASE_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data"
SALMON_DIR="$BASE_DIR/salmon_quant"
COUNTS_DIR="$BASE_DIR/counts"
TMP_DIR="$COUNTS_DIR/tmp"

mkdir -p "$TMP_DIR"

echo '🔄 Extracting column 4 from quant.sf for each sample...'
for QUANT in "$SALMON_DIR"/*/quant.sf; do
  SAMPLE=$(basename "$(dirname "$QUANT")")
  tail -n +2 "$QUANT" | cut -f4 > "$TMP_DIR/${SAMPLE}.count"
done

echo '🔄 Extracting gene IDs from first quant.sf...'
FIRST_QUANT=$(find "$SALMON_DIR" -name quant.sf | head -n 1)
tail -n +2 "$FIRST_QUANT" | cut -f1 > "$TMP_DIR/geneids.txt"

echo '🔄 Combining gene IDs and counts into table...'
paste "$TMP_DIR/geneids.txt" "$TMP_DIR"/*.count > "$TMP_DIR/tmp.out"

echo '🔄 Adding header row using sorted samples.txt...'
SAMPLES_FILE="$BASE_DIR/preproc/samples.txt"
HEADER=$(cat "$SAMPLES_FILE" | sort | paste -s)
echo -e "gene_id\t$HEADER" > "$COUNTS_DIR/mgptestdata_counts.txt"
cat "$TMP_DIR/tmp.out" >> "$COUNTS_DIR/mgptestdata_counts.txt"

echo '✅ Done! Output written to:'
echo "$COUNTS_DIR/mgptestdata_counts.txt"
