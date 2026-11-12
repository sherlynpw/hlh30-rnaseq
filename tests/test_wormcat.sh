#!/usr/bin/env bash
set -euo pipefail

mkdir -p /tmp/wormcat_test /tmp/wormcat_out

# Create a gene list
cat > /tmp/wormcat_test/genes_de.txt <<'EOF'
WBGene00000001
WBGene00000002
WBGene00000003
WBGene00000004
WBGene00000005
WBGene00000006
WBGene00000007
WBGene00000008
WBGene00000009
WBGene00000010
WBGene00000012
WBGene00000013
WBGene00000015
WBGene00000016
WBGene00000017
WBGene00000018
WBGene00000020
WBGene00000021
WBGene00000022
WBGene00000023
EOF

Rscript bin/wormcat.R --genes /tmp/wormcat_test --outdir /tmp/wormcat_out

echo "=== Output ==="
find /tmp/wormcat_out -type f
