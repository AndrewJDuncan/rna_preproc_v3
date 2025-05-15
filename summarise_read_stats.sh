#!/bin/bash

INTERMED_DIR=~/projects/rna_pipeline/mgp_test_data/intermediary_files
PREPROC_DIR=~/projects/rna_pipeline/mgp_test_data/preproc
OUT_CSV="$PREPROC_DIR/read_summary.csv"

echo -e "\n=========  SUMMARY TABLE (CSV)  ========="
printf "%-30s %-15s %-15s %-12s\n" "Sample" "Initial_Reads" "Final_Reads" "%_Retained"
echo "Sample,Initial_Reads,Final_Reads,Percent_Retained" > "$OUT_CSV"

for FINAL_JSON in "$PREPROC_DIR"/*_stats2.json; do
    SAMPLE=$(basename "$FINAL_JSON" | sed 's/_stats2.json//')

    # Try to find matching initial read file
    INITIAL_TXT=$(find "$INTERMED_DIR" -name "${SAMPLE}_initial_read_lines.txt")
    [[ ! -f "$INITIAL_TXT" ]] && INITIAL_TXT=$(find "$INTERMED_DIR" -name "${SAMPLE}_initial_stats.txt")

    if [[ -f "$INITIAL_TXT" ]]; then
        INIT=$(head -n 1 "$INITIAL_TXT" | tr -dc '0-9')
    else
        INIT="MISSING"
    fi

    FINAL=$(jq '.in1_read_count' "$FINAL_JSON" 2>/dev/null)
    if [[ "$INIT" =~ ^[0-9]+$ ]] && [[ "$FINAL" =~ ^[0-9]+$ ]]; then
        RETAINED=$(awk "BEGIN { printf \"%.1f\", ($FINAL/$INIT)*100 }")
    else
        RETAINED="N/A"
    fi

    printf "%-30s %-15s %-15s %-12s\n" "$SAMPLE" "$INIT" "$FINAL" "$RETAINED"
    echo "$SAMPLE,$INIT,$FINAL,$RETAINED" >> "$OUT_CSV"
done

echo -e "\nCSV summary written to: $OUT_CSV"
