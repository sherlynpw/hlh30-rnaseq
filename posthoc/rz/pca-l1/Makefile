# Fig 9B: L1-only PCA, WT vs hlh-30(tm1978). Same canvas (7.5 x 5.5 in) and
# placed width one third of the page (2.5 in) as 9A; axis text at 10pt:
# FONT_SCALE = 10 / (9.6 * 2.5 / 7.5).
PAPER_FONT_SCALE ?= 3.125

.PHONY: paper fetch clean
.DEFAULT_GOAL := paper

S3_BUCKET ?= s3://celegans-pipeline-rayandrew
fetch:
	mkdir -p cache
	aws s3 cp --only-show-errors $(S3_BUCKET)/results/pca-hlh30-l1/exploration/pca/pca_coordinates.csv cache/pca_coordinates.csv
	aws s3 cp --only-show-errors $(S3_BUCKET)/results/pca-hlh30-l1/exploration/pca/variance_explained.txt cache/variance_explained.txt

paper: fetch
	Rscript pca-l1.R --font-scale $(PAPER_FONT_SCALE) --legend genotype --outdir out

clean:
	rm -rf out
