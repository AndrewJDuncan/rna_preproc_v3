#!/bin/bash
set -euo pipefail

# ===== Activate Conda =====
set +u
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools
set -u

# ===== Parameters =====
THREADS=32
SALMON_INDEX="$HOME/projects/rna_pipeline/references/salmon_index"
INPUT_DIR="$HOME/projects/rna_pipeline/mgp_test_data/preproc"
OUTPUT_DIR="$HOME/projects/rna_pipeline/mgp_test_data/salmon_quant"

mkdir -p "$OUTPUT_DIR"

# ===== Run Salmon for Each Sample =====
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  echo "Salmon running quantification <*)))<"
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
for R1 in "$INPUT_DIR"/*_cleaned_R1.fastq.gz; do
  SAMPLE=$(basename "$R1" | sed 's/_R1_cleaned.fastq.gz//')
  R2="$INPUT_DIR/${SAMPLE}_cleaned_R2.fastq.gz"

  echo -e "\n=============================="
  echo "Quantifying sample: $SAMPLE"
  echo "=============================="

  salmon quant \
    -i "$SALMON_INDEX" \
    -l A \
    -1 "$R1" -2 "$R2" \
    -p "$THREADS" \
    --validateMappings \
    -o "$OUTPUT_DIR/$SAMPLE" \
    --gcBias || { echo "Salmon quant failed for $SAMPLE" >&2; exit 1; }

  echo "✅ Finished $SAMPLE"
done
