#!/bin/bash
set -euo pipefail

# ===== Define Directories =====
RAW_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata"
INTER_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files"
PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"

echo -e "Sample\tRaw_Reads\tAfter_PhiX\tAfter_UMI_Extract\tAfter_Cleaning\tAfter_Alignment\tAfter_Dedup"

for R1_FILE in "$RAW_DIR"/*_R1_001.fastq.gz; do
  SAMPLE=$(basename "$R1_FILE" | sed 's/_R1_001.fastq.gz//')

  # Raw input reads
  RAW_READS=$(zcat "$RAW_DIR/${SAMPLE}_R1_001.fastq.gz" | wc -l)
  RAW_READS=$((RAW_READS / 4))

  # After PhiX removal
  PHIX_READS=$(zcat "$INTER_DIR/${SAMPLE}_R1_nophi.fastq.gz" | wc -l)
  PHIX_READS=$((PHIX_READS / 4))

  # After UMI extraction
  UMI_READS=$(zcat "$INTER_DIR/${SAMPLE}_R1_extracted.fastq.gz" | wc -l)
  UMI_READS=$((UMI_READS / 4))

  # After cleaning
  CLEANED_READS=$(zcat "$INTER_DIR/${SAMPLE}_R1_cleaned.fastq.gz" | wc -l)
  CLEANED_READS=$((CLEANED_READS / 4))

  # After alignment (read pairs)
  if [[ -f "$INTER_DIR/${SAMPLE}_aligned.bam" ]]; then
    ALIGNED_READS=$(samtools view -c -f 1 "$INTER_DIR/${SAMPLE}_aligned.bam")  # -f 1 = paired
  else
    ALIGNED_READS="NA"
  fi

  # After deduplication
  if [[ -f "$PREPROC_DIR/${SAMPLE}_dedup.bam" ]]; then
    DEDUP_READS=$(samtools view -c -f 1 "$PREPROC_DIR/${SAMPLE}_dedup.bam")
  else
    DEDUP_READS="NA"
  fi

  echo -e "${SAMPLE}\t${RAW_READS}\t${PHIX_READS}\t${UMI_READS}\t${CLEANED_READS}\t${ALIGNED_READS}\t${DEDUP_READS}"
done
