#!/bin/bash

# ===== Activate Conda Environment =====
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rna-tools || { echo "Failed to activate conda env 'rna-tools'" >&2; exit 1; }

# ===== Config =====
THREADS=32
RAW_DIR=/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/rawdata
INTERMEDIATE_DIR=/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/intermediary_files
PREPROC_DIR=/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc
RRNA_REFERENCE=/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/references/human_rrna.fasta  # <-- Update if needed

mkdir -p "$INTERMEDIATE_DIR" "$PREPROC_DIR"

# ===== Process each sample =====
for R1 in "$RAW_DIR"/*_R1_001.fastq.gz; do
  SAMPLE=$(basename "$R1" | sed 's/_R1_001.fastq.gz//')
  R2="$RAW_DIR/${SAMPLE}_R2_001.fastq.gz"
  JSON_STATS="$INTERMEDIATE_DIR/${SAMPLE}_htsStats.json"
  OUTPUT_PREFIX="$PREPROC_DIR/${SAMPLE}_clean"

  echo -e "\n=============================="
  echo "Processing sample: $SAMPLE"
  echo "=============================="

  # Ensure both R1 and R2 exist
  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "ERROR: Missing input files for $SAMPLE" >&2
    continue
  fi

  echo "[1/10] Running initial stats..."
  hts_Stats -t "$THREADS" -L "$JSON_STATS" -N "initial stats" -1 "$R1" -2 "$R2" -F | \
  tee >(cat >/dev/null) |

  echo "[2/10] Removing PhiX..."
  hts_SeqScreener -A "$JSON_STATS" -N "screen phix" -F | \
  tee >(cat >/dev/null) |

  echo "[3/10] Counting rRNA (not removing)..."
  hts_SeqScreener -A "$JSON_STATS" -N "count rRNA" -r -s "$RRNA_REFERENCE" -F | \
  tee >(cat >/dev/null) |

  echo "[4/10] Removing PCR duplicates..."
  hts_SuperDeduper -A "$JSON_STATS" -N "remove duplicates" -F | \
  tee >(cat >/dev/null) |

  echo "[5/10] Trimming adapters..."
  hts_AdapterTrimmer -A "$JSON_STATS" -N "trim adapters" -F | \
  tee >(cat >/dev/null) |

  echo "[6/10] Trimming polyA/T..."
  hts_PolyATTrim -A "$JSON_STATS" -N "trim polyA/T" -F | \
  tee >(cat >/dev/null) |

  echo "[7/10] Removing Ns..."
  hts_NTrimmer -A "$JSON_STATS" -N "trim Ns" -F | \
  tee >(cat >/dev/null) |

  echo "[8/10] Quality trimming..."
  hts_QWindowTrim -A "$JSON_STATS" -N "Q trimming" -F | \
  tee >(cat >/dev/null) |

  echo "[9/10] Filtering short reads..."
  hts_LengthFilter -A "$JSON_STATS" -N "length filter" -n -m 50 -F | \
  tee >(cat >/dev/null) |

  echo "[10/10] Final stats + write output..."
  hts_Stats -A "$JSON_STATS" -N "final stats" -f "$OUTPUT_PREFIX" -F

  if [[ -f "${OUTPUT_PREFIX}_R1.fastq.gz" && -f "${OUTPUT_PREFIX}_R2.fastq.gz" ]]; then
    echo "✅ Sample $SAMPLE completed successfully."
  else
    echo "❌ ERROR: Output files not created for $SAMPLE" >&2
  fi
done
