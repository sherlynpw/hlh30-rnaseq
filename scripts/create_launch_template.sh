#!/usr/bin/env bash
set -euo pipefail

# Create an EC2 launch template that installs aws CLI on Batch instances
# and makes it available inside containers via bind mount

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

TEMPLATE_NAME="nf-batch-awscli"

# User data script: runs on instance boot, installs aws CLI v2
USER_DATA=$(cat <<'SCRIPT' | base64 -w 0
MIME-Version: 1.0
Content-Type: multipart/mixed; boundary="==BOUNDARY=="

--==BOUNDARY==
Content-Type: text/x-shellscript; charset="us-ascii"

#!/bin/bash
curl -fsSL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "/tmp/awscliv2.zip"
unzip -q /tmp/awscliv2.zip -d /tmp
/tmp/aws/install -i /usr/local/aws-cli -b /usr/local/bin
rm -rf /tmp/awscliv2.zip /tmp/aws

--==BOUNDARY==--
SCRIPT
)

echo "Creating launch template: ${TEMPLATE_NAME}"

aws ec2 create-launch-template \
    --launch-template-name "${TEMPLATE_NAME}" \
    --launch-template-data "{\"UserData\": \"${USER_DATA}\"}" \
    --region "${AWS_REGION}" 2>/dev/null || \
aws ec2 create-launch-template-version \
    --launch-template-name "${TEMPLATE_NAME}" \
    --launch-template-data "{\"UserData\": \"${USER_DATA}\"}" \
    --region "${AWS_REGION}"

echo "Done. Now update the compute environment to use this template."
echo ""
echo "Run:"
echo "  aws batch update-compute-environment --compute-environment ${COMPUTE_ENV} \\"
echo "    --compute-resources '{\"launchTemplate\":{\"launchTemplateName\":\"${TEMPLATE_NAME}\"}}'"
