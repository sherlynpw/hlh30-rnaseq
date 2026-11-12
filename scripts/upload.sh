#!/usr/bin/env bash
set -euo pipefail

# Upload BAM files referenced in samplesheet to S3 (parallel)
# Reads bucket from env.yaml
#
# Usage:
#   ./scripts/upload.sh samples/VZRJXQ_samplesheet.yaml
#   ./scripts/upload.sh samples/VZRJXQ_samplesheet.yaml --dryrun

SAMPLESHEET="${1:?Usage: ./scripts/upload.sh <samplesheet.yaml> [--dryrun]}"
DRYRUN="${2:-}"
MAX_PARALLEL="${MAX_PARALLEL:-4}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

BUCKET="${AWS_BUCKET}"

echo "=== Upload BAMs to s3://${BUCKET} ==="
echo "Samplesheet:  ${SAMPLESHEET}"
echo "Parallel:     ${MAX_PARALLEL} uploads"
[[ -n "$DRYRUN" ]] && echo "Mode:         DRY RUN"
echo ""

upload_one() {
    local s3_path="$1"
    local bucket="$2"
    local dryrun="$3"
    local local_path
    local_path=$(echo "$s3_path" | sed "s|s3://${bucket}/||")

    if [[ ! -f "$local_path" && ! -L "$local_path" ]]; then
        echo "SKIP (not found): $local_path"
        return
    fi

    echo "Uploading: $local_path"
    aws s3 cp "$local_path" "$s3_path" --follow-symlinks $dryrun --quiet

    if [[ -f "${local_path}.bai" || -L "${local_path}.bai" ]]; then
        aws s3 cp "${local_path}.bai" "${s3_path}.bai" --follow-symlinks $dryrun --quiet
    fi
    echo "Done:      $local_path"
}

export -f upload_one

grep 'bam:' "$SAMPLESHEET" | grep -v '^\s*#' | sed 's/.*bam:\s*//' | \
    xargs -P "$MAX_PARALLEL" -I{} bash -c "upload_one '{}' '${BUCKET}' '${DRYRUN}'"

echo ""
echo "=== All uploads complete ==="
