#!/bin/bash

# ===== Activate Conda Environment =====
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools || { echo "Failed to activate conda env 'rna-tools'" >&2; exit 1; }

# === Ensure HTStream tools are available ===
export PATH=~/miniforge3/envs/rna-tools/bin:$PATH

# ===== Config =====
THREADS=32
RAW_DIR=/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata
INTERMEDIATE_DIR=/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files
PREPROC_DIR=/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc

mkdir -p "$PREPROC_DIR" "$INTERMEDIATE_DIR"

for R1 in "$RAW_DIR"/*_R1_001.fastq.gz; do
  SAMPLE=$(basename "$R1" | sed 's/_R1_001.fastq.gz//')
  R2="$RAW_DIR/${SAMPLE}_R2_001.fastq.gz"

  echo -e "\n========="
  echo "Processing sample: $SAMPLE"
  echo "========="

  STATS1="$INTERMEDIATE_DIR/${SAMPLE}_stats1.json"
  STATS2="$INTERMEDIATE_DIR/${SAMPLE}_stats2.json"

  # ========== QC STATS BEFORE ========== 
  hts_Stats -t "$THREADS" -1 "$R1" -2 "$R2" -F > "$STATS1"

  # ========== PHIX REMOVAL ==========
  echo "Running hts_SeqScreener on $SAMPLE"

  PHIX_R1="$INTERMEDIATE_DIR/$(basename "$R1" .fastq.gz)_phix_R1.fastq.gz"
  PHIX_R2="$INTERMEDIATE_DIR/$(basename "$R2" .fastq.gz)_phix_R2.fastq.gz"

  hts_SeqScreener phix \
    -1 "$R1" \
    -2 "$R2" \
    -t "$THREADS" \
    --out1 "$PHIX_R1" \
    --out2 "$PHIX_R2" \
    -F

  if [[ $? -ne 0 ]]; then
    echo "ERROR: hts_SeqScreener failed on $SAMPLE" >&2
    exit 1
  fi

  # ===== Locate generated files =====
  PHIX_R1="$INTERMEDIATE_DIR/$(basename "$R1" .fastq.gz)_phix_R1.fastq.gz"
  PHIX_R2="$INTERMEDIATE_DIR/$(basename "$R2" .fastq.gz)_phix_R2.fastq.gz"

  if [[ ! -f "$PHIX_R1" || ! -f "$PHIX_R2" ]]; then
    echo "ERROR: Expected phix-cleaned files not found for $SAMPLE" >&2
    exit 1
  fi

  # ========== QC STATS AFTER PHIX ==========
  hts_Stats -t "$THREADS" -1 "$PHIX_R1" -2 "$PHIX_R2" -o "$STATS2" -F

  echo "Sample $SAMPLE complete."
done
