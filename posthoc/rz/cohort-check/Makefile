# Room-temperature vs 20 C adult cohorts of the 28-sample set: does the
# cohort change the hlh-30 vs WT DEGs? See cohort-check.R.
#
#   make          # fetch the pipeline inputs, then run in the pipeline container
#   make clean

S3_BUCKET ?= s3://celegans-pipeline-rayandrew
RUN       := VZRJXQ-ruv-rbe
IMG       ?= rayandrew/celegans-r-analysis:0.11.0
R         := $(S3_BUCKET)/results/$(RUN)

.PHONY: all fetch clean
all: out/cohort_summary.tsv

fetch: cache/gene_counts.tsv cache/sample_meta.tsv cache/w_factors.tsv
cache/gene_counts.tsv:
	mkdir -p cache; aws s3 cp --only-show-errors $(R)/counts/gene_counts.tsv $@
cache/sample_meta.tsv:
	mkdir -p cache; aws s3 cp --only-show-errors $(R)/counts/sample_meta.tsv $@
cache/w_factors.tsv:
	mkdir -p cache; aws s3 cp --only-show-errors $(R)/ruvseq/w_factors.tsv $@

out/cohort_summary.tsv: cohort-check.R fetch
	docker run --rm --user $$(id -u):$$(id -g) -v $(abspath ..):/rz -w /rz/cohort-check $(IMG) \
	  Rscript cohort-check.R --cache cache --outdir out --palette ../_common/palette.R

clean:
	rm -rf cache out
