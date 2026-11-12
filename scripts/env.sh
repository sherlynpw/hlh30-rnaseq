#!/usr/bin/env bash
# Source this to load env.yaml values as shell variables
# Usage: source scripts/env.sh

ENV_FILE="${ENV_FILE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/env.yaml}"

if [[ ! -f "$ENV_FILE" ]]; then
    echo "Error: env.yaml not found at $ENV_FILE" >&2
    return 1 2>/dev/null || exit 1
fi

export AWS_ACCOUNT_ID=$(grep 'account_id' "$ENV_FILE" | head -1 | sed 's/.*: *"\?\([^"]*\)"\?.*/\1/')
if [[ -z "$AWS_ACCOUNT_ID" ]]; then
    AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text 2>/dev/null || true)
    export AWS_ACCOUNT_ID
fi
export AWS_REGION=$(grep 'region' "$ENV_FILE" | head -1 | sed 's/.*: *"\?\([^"]*\)"\?/\1/')
export AWS_BUCKET=$(grep 'bucket' "$ENV_FILE" | head -1 | sed 's/.*: *"\?\([^"]*\)"\?/\1/')
export BATCH_QUEUE=$(grep 'batch_queue' "$ENV_FILE" | head -1 | sed 's/.*: *"\?\([^"]*\)"\?/\1/')
export COMPUTE_ENV=$(grep 'compute_env' "$ENV_FILE" | head -1 | sed 's/.*: *"\?\([^"]*\)"\?/\1/')
export MAX_VCPUS=$(grep 'max_vcpus' "$ENV_FILE" | head -1 | sed 's/.*: *"\?\([^"]*\)"\?/\1/')
export ENSEMBL_RELEASE=$(grep 'ensembl_release' "$ENV_FILE" | head -1 | sed 's/.*: *"\?\([^"]*\)"\?/\1/')
export GENOME=$(grep 'genome' "$ENV_FILE" | head -1 | sed 's/.*: *"\?\([^"]*\)"\?/\1/')
