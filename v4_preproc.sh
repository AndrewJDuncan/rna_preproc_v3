#!/bin/bash
set -euo pipefail

# Safely activate Conda
set +u
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools
set -u

# ===== Define Directories =====
RAW_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata"
INTER_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files"
PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"
PHIX_REF="/raid/VIDRL-USERS/HOME/aduncan/bbmap/resources/phix174_ill.ref.fa.gz"

mkdir -p "$INTER_DIR" "$PREPROC_DIR"

# ===== Main Processing Loop =====
for R1_FILE in "$RAW_DIR"/*_R1_001.fastq.gz; do
  SAMPLE=$(basename "$R1_FILE" | sed 's/_R1_001.fastq.gz//')
  R2_FILE="$RAW_DIR/${SAMPLE}_R2_001.fastq.gz"

  echo -e "\n=============================="
  echo "Processing sample: $SAMPLE"
  echo "=============================="

  # ----- Step 1: Initial stats -----
  echo "[1/6] Calculating initial stats..."
  zcat "$R1_FILE" "$R2_FILE" | wc -l > "$INTER_DIR/${SAMPLE}_initial_read_lines.txt"
  echo "$SAMPLE: Raw read count = $(($(cat "$INTER_DIR/${SAMPLE}_initial_read_lines.txt") / 4))" > "$INTER_DIR/${SAMPLE}_initial_stats.txt"

  # ----- Step 2: Remove PhiX using BBduk -----
  echo "[2/6] Removing PhiX..."
  bbduk.sh \
    in1="$R1_FILE" \
    in2="$R2_FILE" \
    out1="$INTER_DIR/${SAMPLE}_R1_nophi.fastq.gz" \
    out2="$INTER_DIR/${SAMPLE}_R2_nophi.fastq.gz" \
    ref="$PHIX_REF" \
    k=31 hdist=1 \
    stats="$INTER_DIR/${SAMPLE}_phix_removal.txt" || { echo "PhiX removal failed for $SAMPLE" >&2; exit 1; }

  # ----- Step 3: Extract UMIs -----
  echo "[3/6] Extracting UMIs..."
  umi_tools extract \
    --extract-method=string \
    --bc-pattern=NNNNNNNN \
    --stdin="$INTER_DIR/${SAMPLE}_R1_nophi.fastq.gz" \
    --stdout="$INTER_DIR/${SAMPLE}_R1_extracted.fastq.gz" \
    --read2-in="$INTER_DIR/${SAMPLE}_R2_nophi.fastq.gz" \
    --read2-out="$INTER_DIR/${SAMPLE}_R2_extracted.fastq.gz" \
    --log="$INTER_DIR/${SAMPLE}_extract.log" || { echo "UMI extraction failed for $SAMPLE" >&2; exit 1; }

  # ----- Step 3.5: Clean FASTQs -----
  echo "[3.5] Cleaning reads (polyA/T trimming, N removal, length filter)..."

  # Trim polyA/T tails (≥10bp A/T at ends), remove Ns, and filter short reads
  bbduk.sh \
    in1="$INTER_DIR/${SAMPLE}_R1_extracted.fastq.gz" \
    in2="$INTER_DIR/${SAMPLE}_R2_extracted.fastq.gz" \
    out1="$INTER_DIR/${SAMPLE}_R1_cleaned.fastq.gz" \
    out2="$INTER_DIR/${SAMPLE}_R2_cleaned.fastq.gz" \
    ref="$REFERENCE_DIR/polyA.fa.gz" \
    k=13 ktrim=r \
    hdist=1 \
    minlength=50 \
    maxns=20 \
    qtrim=r trimq=10 \
    stats="$INTER_DIR/${SAMPLE}_bbduk_cleaning.txt"

  # ----- Step 4: Align to genome (example placeholder) -----
  echo "[4/6] Aligning with HISAT2..."
  hisat2 -p 8 -x /raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references/hg38_index \
  -1 "$INTER_DIR/${SAMPLE}_R1_cleaned.fastq.gz" \
  -2 "$INTER_DIR/${SAMPLE}_R2_cleaned.fastq.gz" \
  | samtools view -bS - > "$INTER_DIR/${SAMPLE}_aligned.bam"

  # Added: sort and index the BAM
  samtools sort -o "$INTER_DIR/${SAMPLE}_aligned_sorted.bam" "$INTER_DIR/${SAMPLE}_aligned.bam"
  mv "$INTER_DIR/${SAMPLE}_aligned_sorted.bam" "$INTER_DIR/${SAMPLE}_aligned.bam"
  samtools index "$INTER_DIR/${SAMPLE}_aligned.bam"

  # ----- Step 5: Deduplicate -----
  echo "[5/6] Deduplicating UMIs..."
  umi_tools dedup \
    --paired \
    -I "$INTER_DIR/${SAMPLE}_aligned.bam" \
    -S "$PREPROC_DIR/${SAMPLE}_dedup.bam" \
    --log="$INTER_DIR/${SAMPLE}_dedup.log" || { echo "UMI deduplication failed for $SAMPLE" >&2; exit 1; }

  # ----- Step 6: Final stats -----
  echo "[6/6] Calculating final stats..."
  samtools flagstat "$PREPROC_DIR/${SAMPLE}_dedup.bam" > "$PREPROC_DIR/${SAMPLE}_flagstat.txt"

  echo "✅ Sample $SAMPLE processed successfully."
done
