#!/bin/bash

# ===== CONFIG: adjust if needed =====
PREPROC_DIR="/raid/VIDRL-USERS/HOME/aduncan/projects/rna_pipeline/mgp_test_data/preproc"

echo ""
echo "========= SUMMARY TABLE ========="
printf "%-30s\t%-15s\t%-15s\t%-10s\n" "Sample" "Initial_Reads" "Final_Reads" "%_Retained"
echo "--------------------------------------------------------------------------"

for stats1 in "$PREPROC_DIR"/*_stats1.json; do
    sample=$(basename "$stats1" | sed 's/_stats1.json//')
    stats2="$PREPROC_DIR/${sample}_stats2.json"

    if [[ -f "$stats2" ]]; then
        initial=$(jq '.total_reads' "$stats1")
        final=$(jq '.total_reads' "$stats2")
        if [[ "$initial" -gt 0 ]]; then
            percent=$(awk "BEGIN {printf \"%.1f\", 100 * $final / $initial}")
        else
            percent="0.0"
        fi
        printf "%-30s\t%-15s\t%-15s\t%-10s\n" "$sample" "$initial" "$final" "$percent"
    else
        printf "%-30s\t%-15s\t%-15s\t%-10s\n" "$sample" "$initial" "MISSING" "N/A"
    fi
done
