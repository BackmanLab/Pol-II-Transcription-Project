#!/bin/bash
# Author: Lucas Carter
# PI: Vadim Backman
# Description: Run after all 00_count_reads.sh array tasks complete.
# Collects per-task count files, computes fractions, writes sample_sheet_fractions.txt.
#
# Usage: bash 00b_make_samplesheet.sh

COUNTS_DIR=counts

# Check all 6 count files exist
echo "Checking count files..."
MISSING=0
for COND in WT_HCT116_CTRL Pol2_6hrs_Aux 1hr_ActD; do
    for REP in 1 2; do
        F=${COUNTS_DIR}/count_${COND}_Rep${REP}.txt
        if [[ ! -f "$F" ]]; then
            echo "  MISSING: $F" >&2
            MISSING=$((MISSING+1))
        else
            echo "  OK: $(cat $F)"
        fi
    done
done

if [[ $MISSING -gt 0 ]]; then
    echo "ERROR: $MISSING count file(s) missing. Wait for all array tasks to finish." >&2
    exit 1
fi

# Concatenate into one file and compute fractions
cat ${COUNTS_DIR}/count_*.txt > sample_sheet.txt

python3 - << 'PYEOF'
import os

rows = []
with open("sample_sheet.txt") as f:
    for line in f:
        parts = line.strip().split('\t')
        cond, rep, r1, r2, n = parts[0], parts[1], parts[2], parts[3], int(parts[4])
        rows.append({'cond': cond, 'rep': rep, 'r1': r1, 'r2': r2, 'n': n})

# Sum per condition
cond_totals = {}
for r in rows:
    cond_totals[r['cond']] = cond_totals.get(r['cond'], 0) + r['n']

min_total = min(cond_totals.values())
min_cond  = min(cond_totals, key=lambda x: cond_totals[x])

print(f"\nRead counts per condition:")
for c, n in sorted(cond_totals.items()):
    print(f"  {c}: {n:,}")

print(f"\nTarget depth: {min_total:,} ({min_cond})")
print()

with open("sample_sheet_fractions.txt", "w") as out:
    out.write("condition\trep\tr1_file\tr2_file\tn_reads\tcond_total\tfraction\n")
    for r in rows:
        frac = min(min_total / cond_totals[r['cond']], 1.0)
        out.write(f"{r['cond']}\t{r['rep']}\t{r['r1']}\t{r['r2']}\t"
                  f"{r['n']}\t{cond_totals[r['cond']]}\t{frac:.6f}\n")
        print(f"  {r['cond']} Rep{r['rep']}: {r['n']:,} reads  fraction={frac:.4f}")

print("\nWritten: sample_sheet_fractions.txt")
PYEOF
