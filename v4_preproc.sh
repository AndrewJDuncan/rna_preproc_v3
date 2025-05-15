#!/bin/bash
set -euo pipefail

# ===== Activate Conda Environment =====
set +u
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools
set -u

# ===== Define Directories =====
RAW_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata"
INTER_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files"
PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"
PHIX_REF="/raid/VIDRL-USERS/HOME/aduncan/bbmap/resources/phix174_ill.ref.fa.gz"
REFERENCE_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references"
POLYA_REF="/raid/VIDRL-USERS/HOME/aduncan/bbmap/resources/polyA.fa.gz"
THREADS=32

mkdir -p "$INTER_DIR" "$PREPROC_DIR"

# ===== Main Processing Loop =====
for R1_FILE in "$RAW_DIR"/*_R1_001.fastq.gz; do
  SAMPLE=$(basename "$R1_FILE" | sed 's/_R1_001.fastq.gz//')
  R2_FILE="$RAW_DIR/${SAMPLE}_R2_001.fastq.gz"

  echo -e "\n=============================="
  echo "Processing sample: $SAMPLE"
  echo "=============================="

  # ----- Step 1: Initial stats -----
  echo "[1/9] Calculating initial stats..."
  zcat "$R1_FILE" "$R2_FILE" | wc -l > "$INTER_DIR/${SAMPLE}_initial_read_lines.txt"
  echo "$SAMPLE: Raw read count = $(($(cat "$INTER_DIR/${SAMPLE}_initial_read_lines.txt") / 4))" > "$INTER_DIR/${SAMPLE}_initial_stats.txt"

  # ----- Step 2: Remove PhiX using BBduk -----
  echo "[2/9] Removing PhiX..."
  bbduk.sh \
    in1="$R1_FILE" \
    in2="$R2_FILE" \
    out1="$INTER_DIR/${SAMPLE}_R1_nophi.fastq.gz" \
    out2="$INTER_DIR/${SAMPLE}_R2_nophi.fastq.gz" \
    ref="$PHIX_REF" \
    k=31 hdist=1 threads=$THREADS \
    stats="$INTER_DIR/${SAMPLE}_phix_removal.txt" || { echo "PhiX removal failed for $SAMPLE" >&2; exit 1; }

  # ----- Step 3: Extract UMIs -----
  echo "[3/9] Extracting UMIs..."
  umi_tools extract \
    --extract-method=string \
    --bc-pattern=NNNNNNNN \
    --stdin="$INTER_DIR/${SAMPLE}_R1_nophi.fastq.gz" \
    --stdout="$INTER_DIR/${SAMPLE}_R1_extracted.fastq.gz" \
    --read2-in="$INTER_DIR/${SAMPLE}_R2_nophi.fastq.gz" \
    --read2-out="$INTER_DIR/${SAMPLE}_R2_extracted.fastq.gz" \
    --log="$INTER_DIR/${SAMPLE}_extract.log" || { echo "UMI extraction failed for $SAMPLE" >&2; exit 1; }
 
   # ----- Step 4: Cleaning reads (adapters, polyA/T, Ns, short reads, quality) -----
   echo "[4/9] Cleaning reads (adapters, polyA/T, N filtering, quality trimming)..."
   bbduk.sh \
      in1="$INTER_DIR/${SAMPLE}_R1_extracted.fastq.gz" \
      in2="$INTER_DIR/${SAMPLE}_R2_extracted.fastq.gz" \
      out1="$INTER_DIR/${SAMPLE}_R1_cleaned.fastq.gz" \
      out2="$INTER_DIR/${SAMPLE}_R2_cleaned.fastq.gz" \
      ref=adapters,"$POLYA_REF" \
      k=23 mink=11 ktrim=r hdist=1 \
      tbo tpe \
      minlength=50 maxns=20 \
      qtrim=r trimq=10 \
      threads=$THREADS \
      stats="$INTER_DIR/${SAMPLE}_bbduk_cleaning.txt" || { echo "Read cleaning failed for $SAMPLE" >&2; exit 1; }
  
   # ----- Step 5: Align to genome -----
   echo "[5/9] Aligning with HISAT2..."
   hisat2 -p 8 -x /raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references/hg38_index \
   -1 "$INTER_DIR/${SAMPLE}_R1_cleaned.fastq.gz" \
   -2 "$INTER_DIR/${SAMPLE}_R2_cleaned.fastq.gz" \
   | samtools view -@ $THREADS -bS - > "$INTER_DIR/${SAMPLE}_aligned.bam"

   # ----- Step 6: sort and index the BAM -----
   echo "[6/9] Sorting and indexing BAM..."
   samtools sort -@ $THREADS -o "$INTER_DIR/${SAMPLE}_aligned_sorted.bam" "$INTER_DIR/${SAMPLE}_aligned.bam"
   mv "$INTER_DIR/${SAMPLE}_aligned_sorted.bam" "$INTER_DIR/${SAMPLE}_aligned.bam"
   samtools index "$INTER_DIR/${SAMPLE}_aligned.bam"

   # ----- Step 7: Deduplicate -----
   echo "[7/9] Deduplicating UMIs..."
   umi_tools dedup \
     --paired \
     -I "$INTER_DIR/${SAMPLE}_aligned.bam" \
     -S "$PREPROC_DIR/${SAMPLE}_dedup.bam" \
     --log="$INTER_DIR/${SAMPLE}_dedup.log" || { echo "UMI deduplication failed for $SAMPLE" >&2; exit 1; }

   # ----- Step 8: Convert deduplicated BAM to FASTQ -----
   echo "[8/9] Converting deduplicated BAM to FASTQ..."
   samtools sort -n -@ $THREADS -o "$INTER_DIR/${SAMPLE}_dedup_name_sorted.bam" "$PREPROC_DIR/${SAMPLE}_dedup.bam"

   bedtools bamtofastq \
     -i "$INTER_DIR/${SAMPLE}_dedup_name_sorted.bam" \
     -fq "$PREPROC_DIR/${SAMPLE}_cleaned_R1.fastq" \
     -fq2 "$PREPROC_DIR/${SAMPLE}_cleaned_R2.fastq" || { echo "FASTQ conversion failed for $SAMPLE" >&2; exit 1; }

   # gzip for Salmon compatibility
   gzip "$PREPROC_DIR/${SAMPLE}_cleaned_R1.fastq"
   gzip "$PREPROC_DIR/${SAMPLE}_cleaned_R2.fastq"
  
   # ----- Step 9: Final stats -----
   echo "[9/9] Final cleaned FASTQ read count..."
   zcat "$INTER_DIR/${SAMPLE}_R1_cleaned.fastq.gz" "$INTER_DIR/${SAMPLE}_R2_cleaned.fastq.gz" | wc -l > "$INTER_DIR/${SAMPLE}_final_read_lines.txt"
   echo "$SAMPLE: Cleaned read count = $(($(cat "$INTER_DIR/${SAMPLE}_final_read_lines.txt") / 4))" > "$INTER_DIR/${SAMPLE}_final_stats.txt"
  

  echo "✅ Sample $SAMPLE processed successfully."
done

 echo -e "\n=============================="
 echo "Generating summary of read retention at each step"
 echo "=============================="

 # Run summary script
 /raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/v2preproc_scripts/summarise_read_stats.sh > "$PREPROC_DIR/read_counts_summary.tsv"

 # Show in terminal:
 column -t "$PREPROC_DIR/read_counts_summary.tsv"
