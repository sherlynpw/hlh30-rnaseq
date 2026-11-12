# HLH-30 integration analyses (Fig S7 to S11 and Supplementary File 6).
# Inputs: ../volcano/cache/hlh30-all-ruv (DESeq2 tables), ../chip-tracks/cache
# (ChIP-Atlas bigWigs, peaks, GTF), published/ (see published/README.md).
#
#   make            # all three
#   make sea        # 1. SEA over JASPAR nematodes + CIS-BP C. elegans on DEG promoters (Docker, ~1 h)
#   make sea-genes  # 1b. the DEGs behind each significant SEA motif (out/sea/sea_motif_genes.xlsx)
#   make bound      # 2+3. HLH-30 binding vs expression; bound genes vs DEGs
#   make published  # 4. our DEGs vs published HLH-30 gene lists; autophagy genes

.PHONY: all fetch sea sea-plot sea-genes bound published
all: sea bound published

# normalized counts and sample sheet of hlh30-all-ruv, for the SEA grid's
# "TF expressed at this stage" check (sea-plot.R)
S3_BUCKET ?= s3://celegans-pipeline-rayandrew
fetch:
	mkdir -p cache
	aws s3 cp --only-show-errors $(S3_BUCKET)/results/hlh30-all-ruv/deseq2/normalized_counts.csv cache/normalized_counts.csv
	aws s3 cp --only-show-errors $(S3_BUCKET)/results/hlh30-all-ruv/counts/sample_meta.tsv cache/sample_meta.tsv

sea:
	Rscript sea-dge.R
	Rscript sea-plot.R
	Rscript sea-genes.R

sea-plot:
	Rscript sea-plot.R

sea-genes:
	Rscript sea-genes.R

bound:
	Rscript bound-vs-dge.R

published:
	Rscript published-overlap.R
