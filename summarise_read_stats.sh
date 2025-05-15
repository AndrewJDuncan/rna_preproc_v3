#!/bin/bash

INTERMED_DIR=~/projects/rna_pipeline/mgp_test_data/intermediary_files
PREPROC_DIR=~/projects/rna_pipeline/mgp_test_data/preproc
OUT_CSV="$PREPROC_DIR/read_summary.csv"

echo -e "\n=========  SUMMARY TABLE (CSV)  ========="
printf "%-30s %-15s %-15s %-12s\n" "Sample" "Initial_Reads" "Final_Reads" "%_Retained"
echo "Sample,Initial_Reads,Final_Reads,Percent_Retained" > "$OUT_CSV"

shopt -s nullglob
for INITIAL_FILE in "$INTERMED_DIR"/*_initial_read_lines.txt; do
    BASENAME=$(basename "$INITIAL_FILE")
    SAMPLE="${BASENAME%_initial_read_lines.txt}"

    INIT=$(head -n 1 "$INITIAL_FILE" | tr -dc '0-9')
    FINAL_FILE="$INTERMED_DIR/${SAMPLE}_final_read_lines.txt"

    if [[ -f "$FINAL_FILE" ]]; then
        FINAL=$(head -n 1 "$FINAL_FILE" | tr -dc '0-9')
    else
        FINAL="MISSING"
    fi

    if [[ "$INIT" =~ ^[0-9]+$ ]] && [[ "$FINAL" =~ ^[0-9]+$ ]]; then
        INIT_READS=$((INIT / 4))
        FINAL_READS=$((FINAL / 4))
        RETAINED=$(awk "BEGIN { printf \"%.1f\", 100 * $FINAL_READS / $INIT_READS }")
    else
        INIT_READS="N/A"
        FINAL_READS="N/A"
        RETAINED="N/A"
    fi

    printf "%-30s %-15s %-15s %-12s\n" "$SAMPLE" "$INIT_READS" "$FINAL_READS" "$RETAINED"
    echo "$SAMPLE,$INIT_READS,$FINAL_READS,$RETAINED" >> "$OUT_CSV"
done
shopt -u nullglob

echo -e "\nCSV summary written to: $OUT_CSV"
