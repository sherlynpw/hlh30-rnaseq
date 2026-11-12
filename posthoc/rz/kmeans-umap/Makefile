# Re-render the pipeline's UMAP panels at paper type size (the pipeline saves
# only the plot, so the embedding is recomputed with the same seed; see the
# script header). One target per panel; condition_from/condition_order are
# copied from each run's samples/<run>.yaml and change the embedding, so they
# are pinned here.
#
#   make            # all four panels into out/fig-*/

# UMAP family: every UMAP panel is placed 4.0 in wide; 10 / (9.6 * 4.0 / 8)
FONT_SCALE ?= 2.083
# Widths follow the 10/03 draft: 3A and 3E half page (3.75 in), 9G one third
# (2.5 in), S1B 4.83 in. Rounded up so 9.6 x scale is a whole point:
FS_3A  ?= 2.2917   # 22 pt in R, 10.3 pt printed at 3.75 in
FS_3E  ?= 2.2917
FS_9G  ?= 3.3333   # 32 pt, 10.0 pt at 2.5 in
FS_S1B ?= 1.7708   # 17 pt, 10.3 pt at 4.83 in
POINT_SIZE ?= 0.6
# 16k genes: smaller still, or the blob is solid colour
POINT_SIZE_9G ?= 0.25

.PHONY: all fig-3a fig-3e fig-9g fig-s1b fetch clean
.DEFAULT_GOAL := all

all: fig-3a fig-3e fig-9g fig-s1b

# inputs from the pipeline results: VST, sample sheet, k-means labels and the
# DEG lists that k-means clustered, per run, plus the HDBSCAN labels for fig
# S1B (cache/VZRJXQ-ruv-rbe-hdb)
S3_BUCKET ?= s3://celegans-pipeline-rayandrew
RUNS := VZRJXQ-ruv-rbe VZRJXQ_no_l4-ruv-rbe hlh30-all-ruv
fetch:
	@for r in $(RUNS); do \
	  mkdir -p cache/$$r; \
	  aws s3 cp --only-show-errors $(S3_BUCKET)/results/$$r/deseq2/vst_counts.csv cache/$$r/vst_counts.csv; \
	  aws s3 cp --only-show-errors $(S3_BUCKET)/results/$$r/counts/sample_meta.tsv cache/$$r/sample_meta.tsv; \
	  aws s3 cp --only-show-errors $(S3_BUCKET)/results/$$r/clustering/kmeans/cluster_assignments.csv cache/$$r/cluster_assignments.csv; \
	  aws s3 cp --only-show-errors --recursive $(S3_BUCKET)/results/$$r/filtered/filtered/ cache/$$r/de-genes/ --exclude '*' --include '*_genes_de.txt'; \
	  Rscript ../_common/renumber-clusters.R --run $$r --csv cache/$$r/cluster_assignments.csv; \
	done
	mkdir -p cache/counts cache/VZRJXQ-ruv-rbe-hdb
	aws s3 cp --only-show-errors $(S3_BUCKET)/results/hlh30-all-ruv/counts/gene_counts.tsv cache/counts/gene_counts.tsv
	aws s3 cp --only-show-errors $(S3_BUCKET)/results/VZRJXQ-ruv-rbe/clustering/hdbscan/cluster_assignments.csv cache/VZRJXQ-ruv-rbe-hdb/cluster_assignments.csv

# fig 3A: VZRJXQ-ruv-rbe, k = 6, placed 3.75 in
fig-3a:
	Rscript kmeans-umap.R --vst cache/VZRJXQ-ruv-rbe/vst_counts.csv \
	  --meta cache/VZRJXQ-ruv-rbe/sample_meta.tsv \
	  --clusters cache/VZRJXQ-ruv-rbe/cluster_assignments.csv \
	  --condition-from genotype_stage \
	  --condition-order "N2_L4,N2_Day3,N2_Day6,hlh-30_L4,hlh-30_Day3,hlh-30_Day6" \
	  --font-scale $(FS_3A) --point-size $(POINT_SIZE) --outdir out/fig-3a

# fig 3E: VZRJXQ_no_l4-ruv-rbe, k = 4, placed 3.75 in
fig-3e:
	Rscript kmeans-umap.R --vst cache/VZRJXQ_no_l4-ruv-rbe/vst_counts.csv \
	  --meta cache/VZRJXQ_no_l4-ruv-rbe/sample_meta.tsv \
	  --clusters cache/VZRJXQ_no_l4-ruv-rbe/cluster_assignments.csv \
	  --condition-from genotype_stage \
	  --condition-order "N2_Day3,N2_Day6,hlh-30_Day3,hlh-30_Day6" \
	  --font-scale $(FS_3E) --point-size $(POINT_SIZE) --outdir out/fig-3e

# fig 9G: hlh30-all-ruv, k = 6 incl. L1, placed 2.5 in; as tall as 9H (3.68 in placed)
fig-9g:
	Rscript kmeans-umap.R --vst cache/hlh30-all-ruv/vst_counts.csv \
	  --meta cache/hlh30-all-ruv/sample_meta.tsv \
	  --clusters cache/hlh30-all-ruv/cluster_assignments.csv \
	  --condition-from genotype_stage \
	  --condition-order "WT_L1,WT_L4,WT_Day3,WT_Day6,tm1978_L1,tm1978_L4,tm1978_Day3,tm1978_Day6" \
	  --title "K-means (k = 6)" --legend bottom --height 11.8 \
	  --font-scale $(FS_9G) --point-size $(POINT_SIZE_9G) --outdir out/fig-9g

# fig S1B: VZRJXQ-ruv-rbe HDBSCAN, placed 4.83 in (point size fixed in script)
fig-s1b:
	Rscript kmeans-umap.R --vst cache/VZRJXQ-ruv-rbe/vst_counts.csv \
	  --meta cache/VZRJXQ-ruv-rbe/sample_meta.tsv \
	  --clusters cache/VZRJXQ-ruv-rbe-hdb/cluster_assignments.csv \
	  --condition-from genotype_stage \
	  --condition-order "N2_L4,N2_Day3,N2_Day6,hlh-30_L4,hlh-30_Day3,hlh-30_Day6" \
	  --mode hdbscan --font-scale $(FS_S1B) --outdir out/fig-s1b

clean:
	rm -rf out
