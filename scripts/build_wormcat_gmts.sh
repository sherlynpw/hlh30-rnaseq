#!/usr/bin/env bash
set -euo pipefail

INPUT_CSV="${1:-whole_genome_v2_nov-11-2021.csv}"
OUTDIR="${2:-.}"
S3_PREFIX="${3:-s3://celegans-pipeline-rayandrew/refs/gene_sets}"

python - "$INPUT_CSV" "$OUTDIR" <<'PY'
import csv
import os
import sys
from collections import defaultdict

input_csv, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

cols = ["Category 1", "Category 2", "Category 3"]
sets = [defaultdict(set) for _ in cols]

with open(input_csv, newline="") as f:
    for row in csv.DictReader(f):
        wb = (row.get("Wormbase ID") or "").strip()
        if not wb:
            continue
        for i, col in enumerate(cols):
            v = (row.get(col) or "").strip()
            if v and v.lower() != "unassigned":
                sets[i][v].add(wb)

for i, col in enumerate(cols):
    name = f"wormcat_cat{i+1}"
    out = os.path.join(outdir, f"{name}.gmt")
    with open(out, "w", newline="") as fh:
        for term in sorted(sets[i]):
            genes = sorted(sets[i][term])
            fh.write("\t".join([term, f"WormCat {col}", *genes]) + "\n")
    print(f"Wrote {len(sets[i])} sets to {out}")
PY

for i in 1 2 3; do
  aws s3 cp "${OUTDIR}/wormcat_cat${i}.gmt" "${S3_PREFIX}/wormcat_cat${i}.gmt"
done
aws s3 ls "${S3_PREFIX}/" | grep wormcat_cat
