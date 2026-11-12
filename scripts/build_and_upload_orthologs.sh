#!/usr/bin/env bash
set -euo pipefail

OUTFILE="${1:-wormbase_to_human_symbol.tsv}"
S3_PATH="${2:-s3://celegans-pipeline-rayandrew/refs/orthologs/wormbase_to_human_symbol.tsv}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TMP_R="$(mktemp "${SCRIPT_DIR}/build_orthologs.XXXXXX.R")"
trap 'rm -f "$TMP_R"' EXIT

cat > "$TMP_R" <<'EOF'
library(biomaRt)

options(timeout = 120)

args <- commandArgs(trailingOnly = TRUE)
outfile <- args[[1]]

mart <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "celegans_gene_ensembl",
  host = "https://useast.ensembl.org"
)

feature_map <- getBM(
  attributes = c("ensembl_gene_id", "wormbase_gene"),
  mart = mart
)

homolog_map <- getBM(
  attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"),
  mart = mart
)

orth <- merge(feature_map, homolog_map, by = "ensembl_gene_id", all = FALSE)
colnames(orth) <- c("ensembl_gene_id", "wormbase_gene", "human_symbol")
orth <- orth[orth$wormbase_gene != "" & orth$human_symbol != "", c("wormbase_gene", "human_symbol"), drop = FALSE]
orth <- unique(orth)
orth <- orth[order(orth$wormbase_gene, orth$human_symbol), , drop = FALSE]

write.table(
  orth,
  outfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("Wrote %d rows to %s\n", nrow(orth), outfile))
EOF

echo "Building ortholog table: ${OUTFILE}"
Rscript "$TMP_R" "$OUTFILE"

echo "Uploading to: ${S3_PATH}"
aws s3 cp "$OUTFILE" "$S3_PATH"

echo "Verifying upload:"
aws s3 ls "$(dirname "$S3_PATH")/"
