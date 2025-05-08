#!/bin/bash

set -euo pipefail

# SAFELY initialize LD_LIBRARY_PATH to avoid 'unbound variable' error
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# Activate the conda environment
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools

# ========== CONFIGURATION ==========
THREADS=32
RAW_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata"
PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"
INTERMEDIATE_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files"

mkdir -p "$PREPROC_DIR" "$INTERMEDIATE_DIR"

# ========== PROCESS EACH SAMPLE ==========
for R1 in "$RAW_DIR"/*_R1_001.fastq.gz; do
    BASENAME=$(basename "$R1" _R1_001.fastq.gz)
    R2="$RAW_DIR/${BASENAME}_R2_001.fastq.gz"
    SAMPLE="$BASENAME"

    echo -e "\n=========\nProcessing sample: $SAMPLE\n========="

    STATS1="$INTERMEDIATE_DIR/${SAMPLE}_stats1.json"
    STATS2="$PREPROC_DIR/${SAMPLE}_stats2.json"
    PHIX_R1="$INTERMEDIATE_DIR/${SAMPLE}_phix_R1.fastq.gz"
    PHIX_R2="$INTERMEDIATE_DIR/${SAMPLE}_phix_R2.fastq.gz"

    # ========== QC STATS BEFORE ==========
    hts_Stats -t "$THREADS" -1 "$R1" -2 "$R2" -F > "$STATS1"

    # ========== PHIX REMOVAL ==========
    hts_SeqScreener phix -1 "$R1" -2 "$R2" -t "$THREADS" "$PHIX_R1" "$PHIX_R2"

    # ========== QC STATS AFTER (PHIX) ==========
    hts_Stats -t "$THREADS" -1 "$PHIX_R1" -2 "$PHIX_R2" -F > "$STATS2"

    # Add more processing steps here as needed, following the same pattern.

    echo "Sample $SAMPLE complete."
done
