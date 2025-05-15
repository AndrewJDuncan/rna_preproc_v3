#!/bin/bash

# ===== CONFIG =====
PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"
CSV_OUTFILE="$PREPROC_DIR/read_summary.csv"

# ===== OUTPUT HEADER =====
echo -e "\n========= SUMMARY TABLE (CSV) ========="
printf "%-30s\t%-15s\t%-15s\t%-10s\n" "Sample" "Initial_Reads" "Final_Reads" "%_Retained"
echo "Sample,Initial_Reads,Final_Reads,Percent_Retained" > "$CSV_OUTFILE"

# ===== MAIN LOOP =====
for stats1 in "$PREPROC_DIR"/*_stats1.json; do
    [[ -e "$stats1" ]] || continue  # skip if no files match
    sample=$(basename "$stats1" _stats1.json)
    stats2="$PREPROC_DIR/${sample}_stats2.json"

    if [[ -f "$stats2" ]]; then
        initial=$(jq -r '.total_reads' "$stats1")
        final=$(jq -r '.total_reads' "$stats2")
        if [[ "$initial" -gt 0 ]]; then
            percent=$(awk "BEGIN {printf \"%.1f\", 100 * $final / $initial}")
        else
            percent="0.0"
        fi
    else
        initial=$(jq -r '.total_reads' "$stats1")
        final="MISSING"
        percent="N/A"
    fi

    printf "%-30s\t%-15s\t%-15s\t%-10s\n" "$sample" "$initial" "$final" "$percent"
    echo "$sample,$initial,$final,$percent" >> "$CSV_OUTFILE"
done

echo -e "\nCSV summary written to: $CSV_OUTFILE"
