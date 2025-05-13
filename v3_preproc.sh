#!/bin/bash

# ===== Activate Conda Environment =====
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools || { echo "Failed to activate conda env 'rna-tools'" >&2; exit 1; }

# ===== Ensure HTStream tools are available =====
export PATH=~/miniforge3/envs/rna-tools/bin:$PATH

# ===== Config =====
THREADS=32
RAW_DIR=/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata
INTERMEDIATE_DIR=/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files
PREPROC_DIR=/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc
RRNA_REFERENCE=/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references/human_rrna.fasta  # <-- Update if different

mkdir -p "$PREPROC_DIR" "$INTERMEDIATE_DIR"

# ===== Main Loop =====
for R1 in "$RAW_DIR"/*_R1_001.fastq.gz; do
  SAMPLE=$(basename "$R1" | sed 's/_R1_001.fastq.gz//')
  R2="$RAW_DIR/${SAMPLE}_R2_001.fastq.gz"

  echo -e "\n========="
  echo "Processing sample: $SAMPLE"
  echo "========="

  STATS1="$INTERMEDIATE_DIR/${SAMPLE}_stats1.json"
  STATS2="$INTERMEDIATE_DIR/${SAMPLE}_stats2.json"
  PHIX_PREFIX="$INTERMEDIATE_DIR/${SAMPLE}_phixclean"

  # ========== Stats + PhiX Removal ==========
  echo "Running hts_Stats + hts_SeqScreener (PhiX removal)"
  hts_Stats -t "$THREADS" -1 "$R1" -2 "$R2" -F -L "$STATS1" |
  hts_SeqScreener phix \
    -A "$STATS1" \
    -f "$PHIX_PREFIX" \
    -F

  PHIX_R1="${PHIX_PREFIX}_R1.fastq.gz"
  PHIX_R2="${PHIX_PREFIX}_R2.fastq.gz"

  if [[ ! -f "$PHIX_R1" || ! -f "$PHIX_R2" ]]; then
    echo "ERROR: Expected phix-cleaned files not found for $SAMPLE" >&2
    exit 1
  fi

  # ========== Optional: rRNA Screening ==========
  echo "Screening for rRNA (no filtering)"
  hts_SeqScreener \
    -r \
    -s "$RRNA_REFERENCE" \
    -f "$INTERMEDIATE_DIR/${SAMPLE}_rrna" \
    -A "$STATS1" \
    -F

  # ========== Stats After PhiX Cleaning ==========
  echo "Running hts_Stats (post-PhiX cleanup)"
  hts_Stats -t "$THREADS" -1 "$PHIX_R1" -2 "$PHIX_R2" -F -L "$STATS2"

  echo "Sample $SAMPLE complete."
done
