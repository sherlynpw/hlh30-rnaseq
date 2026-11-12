#!/usr/bin/env bash
set -euo pipefail

# AWS Batch setup for C. elegans pipeline
# Reads configuration from env.yaml
#
# Usage:
#   ./scripts/setup_aws.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

REGION="${AWS_REGION}"
BUCKET="${AWS_BUCKET}"
ACCOUNT_ID="${AWS_ACCOUNT_ID}"
EMAIL="${EMAIL:-}"

echo "=== C. elegans Pipeline - AWS Setup ==="
echo "Region:     ${REGION}"
echo "Bucket:     ${BUCKET}"
echo "Account:    ${ACCOUNT_ID}"
echo "Max vCPUs:  ${MAX_VCPUS}"
echo ""

# --- Step 1: S3 Bucket ---
echo "--- Step 1: Creating S3 bucket ---"
if aws s3api head-bucket --bucket "${BUCKET}" 2>/dev/null; then
    echo "Bucket already exists: s3://${BUCKET}"
else
    aws s3 mb "s3://${BUCKET}" --region "${REGION}"
    echo "Created: s3://${BUCKET}"
fi

aws s3api put-bucket-lifecycle-configuration \
    --bucket "${BUCKET}" \
    --lifecycle-configuration '{
        "Rules": [{
            "ID": "cleanup-work-dir",
            "Filter": {"Prefix": "work/"},
            "Status": "Enabled",
            "Expiration": {"Days": 7}
        }]
    }'
echo "Lifecycle rule: work/ auto-deletes after 7 days"

# --- Step 2: IAM Roles ---
echo ""
echo "--- Step 2: Creating IAM roles ---"

if ! aws iam get-role --role-name NfBatchInstanceRole 2>/dev/null; then
    aws iam create-role \
        --role-name NfBatchInstanceRole \
        --assume-role-policy-document '{
            "Version": "2012-10-17",
            "Statement": [{
                "Effect": "Allow",
                "Principal": {"Service": "ec2.amazonaws.com"},
                "Action": "sts:AssumeRole"
            }]
        }'
    aws iam attach-role-policy --role-name NfBatchInstanceRole \
        --policy-arn arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role
    aws iam attach-role-policy --role-name NfBatchInstanceRole \
        --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess
    echo "Created: NfBatchInstanceRole"
else
    echo "Exists: NfBatchInstanceRole"
fi

if ! aws iam get-instance-profile --instance-profile-name NfBatchInstanceProfile 2>/dev/null; then
    aws iam create-instance-profile --instance-profile-name NfBatchInstanceProfile
    aws iam add-role-to-instance-profile \
        --instance-profile-name NfBatchInstanceProfile \
        --role-name NfBatchInstanceRole
    echo "Created: NfBatchInstanceProfile"
    echo "Waiting 10s for profile propagation..."
    sleep 10
else
    echo "Exists: NfBatchInstanceProfile"
fi

if ! aws iam get-role --role-name NfBatchServiceRole 2>/dev/null; then
    aws iam create-role \
        --role-name NfBatchServiceRole \
        --assume-role-policy-document '{
            "Version": "2012-10-17",
            "Statement": [{
                "Effect": "Allow",
                "Principal": {"Service": "batch.amazonaws.com"},
                "Action": "sts:AssumeRole"
            }]
        }'
    aws iam attach-role-policy --role-name NfBatchServiceRole \
        --policy-arn arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole
    echo "Created: NfBatchServiceRole"
else
    echo "Exists: NfBatchServiceRole"
fi

if ! aws iam get-role --role-name NfJobRole 2>/dev/null; then
    aws iam create-role \
        --role-name NfJobRole \
        --assume-role-policy-document '{
            "Version": "2012-10-17",
            "Statement": [{
                "Effect": "Allow",
                "Principal": {"Service": "ecs-tasks.amazonaws.com"},
                "Action": "sts:AssumeRole"
            }]
        }'
    POLICY_ARN=$(aws iam create-policy \
        --policy-name NfJobPolicy \
        --policy-document '{
            "Version": "2012-10-17",
            "Statement": [{
                "Effect": "Allow",
                "Action": ["batch:*", "s3:*", "ecs:DescribeTasks",
                           "ec2:DescribeInstances", "logs:*"],
                "Resource": "*"
            }]
        }' \
        --query Policy.Arn --output text 2>/dev/null || \
        echo "arn:aws:iam::${ACCOUNT_ID}:policy/NfJobPolicy")
    aws iam attach-role-policy --role-name NfJobRole --policy-arn "${POLICY_ARN}"
    echo "Created: NfJobRole"
else
    echo "Exists: NfJobRole"
fi

# --- Step 3: AWS Batch ---
echo ""
echo "--- Step 3: Creating Batch compute environment ---"

VPC_ID=$(aws ec2 describe-vpcs --filters Name=isDefault,Values=true \
    --query 'Vpcs[0].VpcId' --output text --region "${REGION}")
SUBNETS=$(aws ec2 describe-subnets --filters Name=vpc-id,Values="${VPC_ID}" \
    --query 'Subnets[*].SubnetId' --output json --region "${REGION}")
SG_ID=$(aws ec2 describe-security-groups \
    --filters Name=vpc-id,Values="${VPC_ID}" Name=group-name,Values=default \
    --query 'SecurityGroups[0].GroupId' --output text --region "${REGION}")

if ! aws batch describe-compute-environments --compute-environments "${COMPUTE_ENV}" \
    --query 'computeEnvironments[0].status' --output text 2>/dev/null | grep -q VALID; then
    aws batch create-compute-environment \
        --compute-environment-name "${COMPUTE_ENV}" \
        --type MANAGED \
        --compute-resources "{
            \"type\": \"SPOT\",
            \"allocationStrategy\": \"SPOT_CAPACITY_OPTIMIZED\",
            \"bidPercentage\": 60,
            \"minvCpus\": 0,
            \"maxvCpus\": ${MAX_VCPUS},
            \"desiredvCpus\": 0,
            \"instanceTypes\": [\"m5\", \"m5a\", \"r5\", \"r5a\", \"c5\", \"c5a\"],
            \"subnets\": ${SUBNETS},
            \"securityGroupIds\": [\"${SG_ID}\"],
            \"instanceRole\": \"arn:aws:iam::${ACCOUNT_ID}:instance-profile/NfBatchInstanceProfile\"
        }" \
        --service-role "arn:aws:iam::${ACCOUNT_ID}:role/NfBatchServiceRole" \
        --region "${REGION}"
    echo "Created: ${COMPUTE_ENV}"
    sleep 30
else
    echo "Exists: ${COMPUTE_ENV}"
fi

if ! aws batch describe-job-queues --job-queues "${BATCH_QUEUE}" \
    --query 'jobQueues[0].status' --output text 2>/dev/null | grep -q VALID; then
    aws batch create-job-queue \
        --job-queue-name "${BATCH_QUEUE}" \
        --priority 1 \
        --compute-environment-order "[{\"order\":1,\"computeEnvironment\":\"${COMPUTE_ENV}\"}]" \
        --region "${REGION}"
    echo "Created: ${BATCH_QUEUE}"
else
    echo "Exists: ${BATCH_QUEUE}"
fi

echo ""
echo "=== Setup complete ==="
