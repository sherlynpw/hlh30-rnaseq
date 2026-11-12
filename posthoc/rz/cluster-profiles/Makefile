# Render k-means cluster profiles (one panel per cluster) from a pipeline run.
#
#   make                  # download + plot for RUN=hlh30-all-ruv
#   make RUN=<other-run>
#   make upload           # sync out/ back to S3
#   make clean | distclean

RUN       ?= hlh30-all-ruv
# Fig 3B/3F placed at ~4.45in of a 9.5in canvas; 9H at 4.34in. See --font-scale.
FONT_SCALE ?= 1
NCOL  ?= 3
WIDTH ?= 9.5
# `make paper`: figs 3B / 3F / 9H at 10pt on the 9.5 in canvas,
# values measured from the rendered PDFs (make measure in ../paper-figures),
# placed 4.41 / 4.48 / 4.34 in.
PAPER_FONT_SCALE_3B ?= 2.555
# 3F is 2 x 2 on a 6.6 in canvas: 10 / (8.8 * 4.48 / 6.6) rounded up
PAPER_FONT_SCALE_3F ?= 1.78
PAPER_FONT_SCALE_9H ?= 2.2438   # 9H two thirds of the page (5.0 in): 2.585 x 4.34 / 5.0
SHORT_STAGES ?=   # set to --short-stages at 10pt
S3_BUCKET ?= s3://celegans-pipeline-rayandrew

S3_VST      := $(S3_BUCKET)/results/$(RUN)/deseq2/vst_counts.csv
S3_META     := $(S3_BUCKET)/results/$(RUN)/counts/sample_meta.tsv
S3_CLUSTERS := $(S3_BUCKET)/results/$(RUN)/clustering/kmeans/cluster_assignments.csv

CACHE  := cache/$(RUN)
VST      := $(CACHE)/vst_counts.csv
META     := $(CACHE)/sample_meta.tsv
CLUSTERS := $(CACHE)/cluster_assignments.csv

LABELS      ?= labels-$(RUN).tsv
PLOT_SCRIPT := cluster-profiles.R
OUTDIR      := out/$(RUN)
STAMP       := $(OUTDIR)/.done

.PHONY: all paper download deps clean distclean print upload
.DEFAULT_GOAL := all

all: $(STAMP)

paper:
	$(MAKE) -B RUN=VZRJXQ-ruv-rbe FONT_SCALE=$(PAPER_FONT_SCALE_3B) SHORT_STAGES=--short-stages
	$(MAKE) -B RUN=VZRJXQ_no_l4-ruv-rbe FONT_SCALE=$(PAPER_FONT_SCALE_3F) SHORT_STAGES=--short-stages NCOL=2 WIDTH=6.6
	$(MAKE) -B RUN=hlh30-all-ruv FONT_SCALE=$(PAPER_FONT_SCALE_9H) SHORT_STAGES=--short-stages

download: $(VST) $(META) $(CLUSTERS)

$(VST):
	mkdir -p $(CACHE)
	aws s3 cp $(S3_VST) $@

$(META):
	mkdir -p $(CACHE)
	aws s3 cp $(S3_META) $@

# pipeline k-means numbers -> paper numbers (../_common/kmeans-renumber.tsv)
$(CLUSTERS):
	mkdir -p $(CACHE)
	aws s3 cp $(S3_CLUSTERS) $@
	Rscript ../_common/renumber-clusters.R --run $(RUN) --csv $@

# LABELS is optional; passed only when the file exists.
$(STAMP): $(VST) $(META) $(CLUSTERS) $(PLOT_SCRIPT)
	mkdir -p $(OUTDIR)
	Rscript $(PLOT_SCRIPT) \
		--vst      $(VST) \
		--meta     $(META) \
		--clusters $(CLUSTERS) \
		$(if $(wildcard $(LABELS)),--labels $(LABELS),) \
		--font-scale $(FONT_SCALE) $(SHORT_STAGES) --ncol $(NCOL) --width $(WIDTH) \
		--outdir   $(OUTDIR)
	touch $@

deps:
	Rscript -e 'install.packages(c("optparse","ggplot2","dplyr","tidyr"), repos="https://cloud.r-project.org")'

S3_UPLOAD_DEST ?= $(S3_BUCKET)/posthoc_results/rz/cluster-profiles

upload:
	aws s3 sync out $(S3_UPLOAD_DEST) --delete

clean:
	rm -rf out

distclean: clean
	rm -rf cache

print:
	@echo RUN=$(RUN)
	@echo LABELS=$(LABELS)
	@echo OUTDIR=$(OUTDIR)
