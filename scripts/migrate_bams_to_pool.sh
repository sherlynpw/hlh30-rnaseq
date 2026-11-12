#!/usr/bin/env bash
set -euo pipefail

# migrate_bams_to_pool.sh - Move per-run nf-core BAM outputs into the shared
# pool layout used by run.sh after the consolidation to approach A.
#
# Before: s3://$AWS_BUCKET/rnaseq_output/<run_id>/{standard,3prime}/star_salmon/*.bam
# After:  s3://$AWS_BUCKET/rnaseq_output/_pool/{standard,3prime}/star_salmon/*.bam
#
# Server-side S3 copies (no download). Existing per-run prefixes are left in
# place; delete them manually once you've verified the pool works for you.
#
# Usage:
#   ./scripts/migrate_bams_to_pool.sh              # dry-run, lists candidates
#   ./scripts/migrate_bams_to_pool.sh --apply      # actually copy
#   ./scripts/migrate_bams_to_pool.sh --apply --delete-source   # also rm originals
#
# Optional: --run-id <id> restricts to a single per-run prefix.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

APPLY=false
DELETE_SOURCE=false
RUN_ID_FILTER=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --apply)          APPLY=true; shift ;;
        --delete-source)  DELETE_SOURCE=true; shift ;;
        --run-id)         RUN_ID_FILTER="$2"; shift 2 ;;
        *)                echo "Unknown option: $1"; exit 1 ;;
    esac
done

POOL="s3://${AWS_BUCKET}/rnaseq_output/_pool"
SRC_ROOT="s3://${AWS_BUCKET}/rnaseq_output"

echo "Source root: ${SRC_ROOT}/"
echo "Pool:        ${POOL}/"
echo "Apply:       ${APPLY}"
echo "Delete src:  ${DELETE_SOURCE}"
[[ -n "${RUN_ID_FILTER}" ]] && echo "Filter:      ${RUN_ID_FILTER}"
echo ""

# List top-level run-id prefixes under SRC_ROOT.
mapfile -t RUN_IDS < <(aws s3 ls "${SRC_ROOT}/" \
    | awk '/PRE /{gsub("/","",$2); print $2}' \
    | grep -v '^_pool$' \
    | { [[ -n "${RUN_ID_FILTER}" ]] && grep -Fx "${RUN_ID_FILTER}" || cat; })

if [[ "${#RUN_IDS[@]}" -eq 0 ]]; then
    echo "No per-run prefixes found to migrate."
    exit 0
fi

for rid in "${RUN_IDS[@]}"; do
    for sub in standard 3prime; do
        src="${SRC_ROOT}/${rid}/${sub}/star_salmon"
        dst="${POOL}/${sub}/star_salmon"
        # Each .markdup.sorted.bam (+ .bai). nf-core publishes both side-by-side.
        mapfile -t bams < <(aws s3 ls "${src}/" 2>/dev/null \
            | awk '{print $NF}' \
            | grep -E '\.markdup\.sorted\.bam(\.bai)?$' || true)
        [[ "${#bams[@]}" -eq 0 ]] && continue
        echo "[${rid}/${sub}] ${#bams[@]} object(s)"
        for b in "${bams[@]}"; do
            if aws s3 ls "${dst}/${b}" >/dev/null 2>&1; then
                echo "  skip (already in pool): ${b}"
                if [[ "${APPLY}" == true && "${DELETE_SOURCE}" == true ]]; then
                    aws s3 rm "${src}/${b}"
                fi
                continue
            fi
            echo "  copy: ${b}"
            if [[ "${APPLY}" == true ]]; then
                aws s3 cp "${src}/${b}" "${dst}/${b}"
                if [[ "${DELETE_SOURCE}" == true ]]; then
                    aws s3 rm "${src}/${b}"
                fi
            fi
        done
    done
done

echo ""
if [[ "${APPLY}" == false ]]; then
    echo "Dry run. Re-run with --apply to actually copy."
fi
