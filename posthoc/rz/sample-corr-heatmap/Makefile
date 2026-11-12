# Sample correlation heatmaps from RUV-corrected VST counts.
#
# For each RUN, two heatmaps are produced:
#   corr_corrected_<METHOD>.pdf  (vst_counts.csv      = post-RUV/RBE)
#   corr_raw_<METHOD>.pdf        (vst_counts_raw.csv  = no batch correction)
#
# Replicates within (genotype, stage) are collapsed by per-gene mean before
# computing the correlation matrix.

S3_BUCKET ?= s3://celegans-pipeline-rayandrew
RUNS      ?= pca-hlh30-all-ruv
METHOD    ?= pearson          # pearson | spearman | kendall
TOP_N     ?= 1000             # 0 = all genes; otherwise top-N most variable

# Canvas pinned to the placed panel size so it goes in 1:1 and 10pt stays 10pt.
# Fig 1C is placed at 5.25 x 4.99 in (133.3 x 126.6 mm), measured with pdfimages.
WIDTH     ?= 0
HEIGHT    ?= 0
FONT_SCALE ?= 1
# `make paper`: the canvas is laid out in points (labels 10pt, cells
# CELL_LINES text lines) so FONT_SCALE 1 is 10pt when placed 1:1. Fig 1C is
# half a page wide (3.75 in, 10/03 draft): 2.22-line cells give exactly
# that width; fig 9C is 1.5-line cells at 1:1 (3.57 in). Both widths are in ../paper-figures/manifest.tsv.
CELL_LINES ?= 1.5
LEGEND     ?= horizontal
MUT_LABEL  ?= hlh-30(tm1978)
CELL_LINES_1C ?= 2.22   # fig 1C half page: 270 pt = 3.75 in at 1:1
FONT_MIN  ?= 10
FIG_SCALE ?= 1

PLOT_SCRIPT := sample-corr-heatmap.R

CACHE  := cache
OUTDIR := out
STAMP  := $(OUTDIR)/.done

.PHONY: all paper deps clean distclean print upload
.DEFAULT_GOAL := all

all: $(STAMP)

paper:
	$(MAKE) -B RUNS=VZRJXQ-ruv-rbe FONT_SCALE=1 CELL_LINES=$(CELL_LINES_1C)
	$(MAKE) -B RUNS=pca-hlh30-all-ruv FONT_SCALE=1

# One stamp covers all runs; rebuild if script or settings change.
$(STAMP): $(PLOT_SCRIPT)
	mkdir -p $(OUTDIR)
	for run in $(RUNS); do \
	  mkdir -p $(CACHE)/$$run $(OUTDIR)/$$run; \
	  for f in vst_counts.csv vst_counts_raw.csv normalized_counts.csv; do \
	    [ -f $(CACHE)/$$run/$$f ] || aws s3 cp $(S3_BUCKET)/results/$$run/deseq2/$$f $(CACHE)/$$run/$$f; \
	  done; \
	  [ -f $(CACHE)/$$run/sample_meta.tsv ] || aws s3 cp $(S3_BUCKET)/results/$$run/counts/sample_meta.tsv $(CACHE)/$$run/sample_meta.tsv; \
	  FONT_MIN=$(FONT_MIN) FIG_SCALE=$(FIG_SCALE) Rscript $(PLOT_SCRIPT) \
	    --vst-corrected  $(CACHE)/$$run/vst_counts.csv \
	    --vst-raw        $(CACHE)/$$run/vst_counts_raw.csv \
	    --log2-norm      $(CACHE)/$$run/normalized_counts.csv \
	    --meta           $(CACHE)/$$run/sample_meta.tsv \
	    --method         $(METHOD) \
	    --top-n          $(TOP_N) \
	    --font-scale     $(FONT_SCALE) \
	    --width          $(WIDTH) \
	    --height         $(HEIGHT) \
	    --cell-lines     $(CELL_LINES) \
	    --legend         $(LEGEND) \
	    --mut-label      "$(MUT_LABEL)" \
	    --outdir         $(OUTDIR)/$$run; \
	done
	touch $@

deps:
	Rscript -e 'install.packages(c("optparse","pheatmap","readr","dplyr","matrixStats"), repos="https://cloud.r-project.org")'

S3_UPLOAD_DEST ?= $(S3_BUCKET)/posthoc_results/rz/sample-corr-heatmap

upload: $(STAMP)
	aws s3 sync $(OUTDIR) $(S3_UPLOAD_DEST) --delete

clean:
	rm -rf $(OUTDIR)

distclean: clean
	rm -rf $(CACHE)

print:
	@echo RUNS=$(RUNS)
	@echo METHOD=$(METHOD)
	@echo TOP_N=$(TOP_N)
	@echo OUTDIR=$(OUTDIR)
