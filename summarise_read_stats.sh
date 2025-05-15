#!/bin/bash
set -euo pipefail

# ===== Activate Conda =====
# Comment out these lines if samtools is already in your PATH
set +u
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools
set -u

# ===== Define Directories =====
RAW_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata"
INTER_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files"
PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"

# ===== Print header =====
printf "%-20s %12s %12s %10s %18s %10s %18s %10s %15s %10s %15s %10s\n" \
"Sample" "Raw_Reads" "After_PhiX" "%_Ret" "After_UMI_Extract" "%_Ret" \
"After_Cleaning" "%_Ret" "Aligned_Reads" "%_Ret" "Dedup_Reads" "%_Ret"

# ===== Loop through samples =====
for R1_FILE in "$RAW_DIR"/*_R1_001.fastq.gz; do
  SAMPLE=$(basename "$R1_FILE" | sed 's/_R1_001.fastq.gz//')

  # ---- Raw reads ----
  RAW_READS=$(zcat "$RAW_DIR/${SAMPLE}_R1_001.fastq.gz" | wc -l)
  RAW_READS=$((RAW_READS / 4))

  # ---- After PhiX ----
  PHIX_READS=$(zcat "$INTER_DIR/${SAMPLE}_R1_nophi.fastq.gz" | wc -l)
  PHIX_READS=$((PHIX_READS / 4))
  PHIX_PCT=$(awk -v a=$RAW_READS -v b=$PHIX_READS 'BEGIN { printf("%.1f", (b/a)*100) }')

  # ---- After UMI extraction ----
  UMI_READS=$(zcat "$INTER_DIR/${SAMPLE}_R1_extracted.fastq.gz" | wc -l)
  UMI_READS=$((UMI_READS / 4))
  UMI_PCT=$(awk -v a=$PHIX_READS -v b=$UMI_READS 'BEGIN { printf("%.1f", (b/a)*100) }')

  # ---- After Cleaning ----
  CLEANED_READS=$(zcat "$INTER_DIR/${SAMPLE}_R1_cleaned.fastq.gz" | wc -l)
  CLEANED_READS=$((CLEANED_READS / 4))
  CLEANED_PCT=$(awk -v a=$UMI_READS -v b=$CLEANED_READS 'BEGIN { printf("%.1f", (b/a)*100) }')

  # ---- After Alignment ----
  if [[ -f "$INTER_DIR/${SAMPLE}_aligned.bam" ]]; then
    ALIGNED_READS=$(samtools view -c -f 1 "$INTER_DIR/${SAMPLE}_aligned.bam")
    TOTAL_CLEANED_READS=$((CLEANED_READS * 2))
    ALIGN_PCT=$(awk -v a=$TOTAL_CLEANED_READS -v b=$ALIGNED_READS 'BEGIN { printf("%.1f", (b/a)*100) }')
  else
    ALIGNED_READS="NA"
    ALIGN_PCT="NA"
  fi

  # ---- After Deduplication ----
  if [[ -f "$PREPROC_DIR/${SAMPLE}_dedup.bam" ]]; then
    DEDUP_READS=$(samtools view -c -f 1 "$PREPROC_DIR/${SAMPLE}_dedup.bam")
    if [[ "$ALIGNED_READS" != "NA" ]]; then
      TOTAL_ALIGNED_READS=$((ALIGNED_READS * 2))
      DEDUP_PCT=$(awk -v a=$TOTAL_ALIGNED_READS -v b=$DEDUP_READS 'BEGIN { printf("%.1f", (b/a)*100) }')
    else
      DEDUP_PCT="NA"
    fi
  else
    DEDUP_READS="NA"
    DEDUP_PCT="NA"
  fi

  # ---- Print result row ----
  printf "%-20s %12s %12s %10s %18s %10s %18s %10s %15s %10s %15s %10s\n" \
  "$SAMPLE" "$RAW_READS" "$PHIX_READS" "${PHIX_PCT}%" "$UMI_READS" "${UMI_PCT}%" \
  "$CLEANED_READS" "${CLEANED_PCT}%" "$ALIGNED_READS" "${ALIGN_PCT}%" "$DEDUP_READS" "${DEDUP_PCT}%"
done
