#!/usr/bin/env bash
# GO with no cutoff on each gene list in lost-terms.tsv, so the terms that are
# no longer significant still get their current (forward-strand) values. Runs
# the pipeline's own bin/go_enrichment.R in its container with the run's own
# universe and gene_db, as the pipeline and compare-de/ do, with
# --pvalue 1 --qvalue 1. Writes out-lost/<run>/<list>/go_<ONT>.csv.
#
#   ./rerun-lost.sh lost-terms.tsv out-lost
set -euo pipefail

lost=$1; out=$2
B=${S3_BUCKET:-s3://celegans-pipeline-rayandrew}
IMG=${IMG:-rayandrew/celegans-r-analysis:0.11.0}
REPO=$(cd "$(dirname "$0")/../../.." && pwd)
work=cache/lost-rerun

grep -v '^#' "$lost" | tail -n +2 | cut -f4,5 | sort -u | while IFS=$'\t' read -r run path; do
  base=$(basename "$path")
  case $path in
    kmeans/*) src=$B/results/$run/clustering/kmeans/$base.txt ;;
    compare_de*/*) src=$B/posthoc_results/rz/compare-de/$run/gene_lists/$base.txt ;;
    *) echo "no gene list source for $path" >&2; exit 1 ;;
  esac
  in=$work/$run/input; mkdir -p "$in" "$out/$run"
  aws s3 cp "$src" "$in/$base.txt" --only-show-errors
  [ -f "$work/$run/background.txt" ] ||
    aws s3 cp "$B/results/$run/enrichment/background.txt" "$work/$run/background.txt" --only-show-errors
done

for run in $(ls "$work"); do
  cp "$REPO/gene_db.sqlite" "$work/$run/gene_db.sqlite"
  docker run --rm --user "$(id -u):$(id -g)" -w /repo -v "$REPO":/repo:ro \
    -v "$(realpath "$work/$run")":/work -v "$(realpath "$out/$run")":/out "$IMG" \
    Rscript bin/go_enrichment.R --input /work/input --gene_db /work/gene_db.sqlite \
    --universe /work/background.txt --pvalue 1 --qvalue 1 --cpus 4 --outdir /out
done
