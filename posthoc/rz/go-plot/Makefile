# Render filtered GO horizontal-bar plots from pipeline enrichment outputs.
#
# Driven by targets.tsv (one row per figure). `make` builds every row.
#   make             # render bar + dot figures for every targets.tsv row
#   make deps        # install required R packages
#   make clean       # remove generated figures
#   make distclean   # also remove cached GO CSVs

S3_BUCKET   ?= s3://celegans-pipeline-rayandrew
CONFIG      ?= targets.tsv
OUTDIR      ?= out
CACHE       ?= cache
# Paper panels are placed at 2.8-4.3 in wide; rendering at the default 7.0 and
# scaling down halves every label. Keep the 7.0 canvas and raise the type:
# `make paper` reads each row's placed_in (measured from the author's layout
# PDF) and scales fonts so the term labels print at TARGET_PT.
WIDTH     ?= 7.0
HEIGHT    ?= 0
WRAP      ?= 38
LEGEND    ?= right
XBREAKS   ?= 5
FONT_SCALE ?= 1
TARGET_PT ?= 0
ROWH      ?= 0.55
FONT_MIN  ?= 10
FIG_SCALE ?= 1
LABEL_FRAC ?= 0

PLOT_SCRIPT := render-bars.R
DOT_SCRIPT  := render-dots.R
COMMON      := go-plot-common.R
STAMP       := $(OUTDIR)/.done

.PHONY: lost all paper full deps clean distclean print
.DEFAULT_GOAL := all

all: $(STAMP)

# The hlh-30 paper panels: one row per figure panel, 10pt in print. Panels
# are half page (3.75 in) or one third (2.5 in, fig 9). Every panel of one
# width has the same label column (0.457 of the width; long terms wrap) and
# the same plot width and grid (0 / step / 2 x step), legend below; only the
# height varies.
paper:
	$(MAKE) CONFIG=targets-paper-figures.tsv OUTDIR=out-paper TARGET_PT=10 ROWH=0.26 \
		LEGEND=bottom LABEL_FRAC=0.457

# Every significant term of each paper panel as fig-<panel>-full.xlsx,
# with the panel's own terms in bold and the terms shown before the
# forward-strand rerun that are no longer significant in grey, with their
# current values from `make lost` (GO with no cutoff on those gene lists).
full: full-tables.R $(COMMON) targets-paper-figures.tsv lost-terms.tsv out-lost/.done
	Rscript full-tables.R --config targets-paper-figures.tsv --outdir out-full \
		--cache $(CACHE) --bucket $(S3_BUCKET) --lost lost-terms.tsv --lost-dir out-lost

lost: out-lost/.done
out-lost/.done: rerun-lost.sh lost-terms.tsv
	S3_BUCKET=$(S3_BUCKET) ./rerun-lost.sh lost-terms.tsv out-lost
	touch $@

$(STAMP): $(CONFIG) $(PLOT_SCRIPT) $(DOT_SCRIPT) $(COMMON)
	mkdir -p $(OUTDIR)
	Rscript $(PLOT_SCRIPT) \
		--config $(CONFIG) \
		--outdir $(OUTDIR) \
		--cache  $(CACHE) \
		--bucket $(S3_BUCKET)
	FONT_MIN=$(FONT_MIN) FIG_SCALE=$(FIG_SCALE) Rscript $(DOT_SCRIPT) \
		--config $(CONFIG) \
		--outdir $(OUTDIR) \
		--cache  $(CACHE) \
		--width  $(WIDTH) \
		--height $(HEIGHT) \
		--font-scale $(FONT_SCALE) \
		--target-pt $(TARGET_PT) \
		--row-height $(ROWH) \
		--wrap   $(WRAP) \
		--legend $(LEGEND) \
		--x-breaks $(XBREAKS) \
		--label-frac $(LABEL_FRAC) \
		--bucket $(S3_BUCKET)
	touch $@

deps:
	Rscript -e 'install.packages(c("optparse","ggplot2","dplyr","readr","stringr","openxlsx2"), repos="https://cloud.r-project.org")'

S3_UPLOAD_DEST ?= $(S3_BUCKET)/posthoc_results/rz/go-plot

upload: $(STAMP)
	aws s3 sync $(OUTDIR) $(S3_UPLOAD_DEST) --delete

clean:
	rm -rf $(OUTDIR) out-full

distclean: clean
	rm -rf $(CACHE)

print:
	@echo CONFIG=$(CONFIG)
	@echo OUTDIR=$(OUTDIR)
	@echo CACHE=$(CACHE)
	@echo S3_BUCKET=$(S3_BUCKET)
