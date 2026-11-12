#!/usr/bin/env bash
set -euo pipefail

# Download Ensembl GTF and upload to S3
# Reads config from env.yaml, updates params file
#
# Usage:
#   ./scripts/upload_refs.sh samples/VZRJXQ-ruv-rbe.yaml

PARAMS="${1:?Usage: ./scripts/upload_refs.sh <params.yaml>}"
TMPDIR="${TMPDIR:-/tmp}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

BUCKET="${AWS_BUCKET}"
RELEASE="${ENSEMBL_RELEASE}"
GTF="Caenorhabditis_elegans.WBcel235.${RELEASE}.gtf"
S3_GTF="s3://${BUCKET}/refs/${GTF}"

echo "=== Upload C. elegans references ==="
echo "Ensembl release: ${RELEASE}"
echo "Destination:     ${S3_GTF}"
echo ""

if [[ ! -f "${TMPDIR}/${GTF}" ]]; then
    echo "Downloading GTF from Ensembl..."
    curl -fSL -o "${TMPDIR}/${GTF}.gz" \
        "https://ftp.ensembl.org/pub/release-${RELEASE}/gtf/caenorhabditis_elegans/${GTF}.gz"
    gunzip -f "${TMPDIR}/${GTF}.gz"
else
    echo "Using cached: ${TMPDIR}/${GTF}"
fi

echo "Uploading to S3..."
aws s3 cp "${TMPDIR}/${GTF}" "${S3_GTF}"

# Update params file
if grep -q 'gtf:.*null' "$PARAMS"; then
    sed -i "s|gtf:.*null|gtf: \"${S3_GTF}\"|" "$PARAMS"
    echo "Updated ${PARAMS} with GTF path."
else
    echo "GTF already set in ${PARAMS}."
fi

echo ""
echo "=== Done ==="
