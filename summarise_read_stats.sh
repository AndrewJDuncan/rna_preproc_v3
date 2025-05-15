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
    FINAL_JSON="$PREPROC_DIR/${SAMPLE}_stats2.json"

    if [[ -f "$FINAL_JSON" ]]; then
        FINAL=$(jq '.in1_read_count' "$FINAL_JSON" 2>/dev/null)
    else
        FINAL="MISSING"
    fi

    if [[ "$INIT" =~ ^[0-9]+$ ]] && [[ "$FINAL" =~ ^[0-9]+$ ]]; then
        RETAINED=$(awk "BEGIN { printf \"%.1f\", ($FINAL/$INIT)*100 }")
    else
        RETAINED="N/A"
    fi

    printf "%-30s %-15s %-15s %-12s\n" "$SAMPLE" "$INIT" "$FINAL" "$RETAINED"
    echo "$SAMPLE,$INIT,$FINAL,$RETAINED" >> "$OUT_CSV"
done
shopt -u nullglob

echo -e "\nCSV summary written to: $OUT_CSV"
