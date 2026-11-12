# Re-render the stage-PCA from a pipeline run's exploration/pca outputs.
#
# Typical use:
#   make                          # download + plot for RUN=VZRJXQ-ruv-rbe
#   make RUN=<other-run>
#   make upload                   # sync out/ back to S3
#   make clean | distclean

RUN       ?= VZRJXQ-ruv-rbe
S3_BUCKET ?= s3://celegans-pipeline-rayandrew

S3_PCA    := $(S3_BUCKET)/results/$(RUN)/exploration/pca
S3_COORDS := $(S3_PCA)/pca_coordinates.csv
S3_VAR    := $(S3_PCA)/variance_explained.txt

# Namespaced by RUN: a bare filename here meant switching RUN silently reused
# the previously downloaded CSV, so out/<RUN>/ could hold another run's plot.
CACHE  := cache/$(RUN)
COORDS := $(CACHE)/pca_coordinates.csv
VAR    := $(CACHE)/variance_explained.txt

# Rendered at the width the panel is *placed* at in the layout, so it goes in
# 1:1 and a 10pt label in R is 10pt in print.
WIDTH     ?= 7.5
HEIGHT    ?= 5.5
# Render at the ORIGINAL canvas (7.5x5.5) and raise FONT_SCALE. Shrinking the
# canvas to the placed size instead doubles every element relative to the plot
# and wrecks the layout -- see fig-1b history.
FONT_SCALE ?= 1
SHORT_STAGES ?=
# `make paper`: fig 1B (placed 3.75 in) and fig 9A (placed 2.77 in), type at
# 10pt: FONT_SCALE = 10 / (9.6 * placed / 7.5).
# Widths follow the 10/03 draft: 1B and 1D half page (3.75 in), 9A and 9B one third (2.5 in)
PAPER_FONT_SCALE_1B ?= 2.0833  # 1B half page, placed 3.75 in: 9.6 x 2.0833 = 20 pt in R, 10.0 pt printed
PAPER_FONT_SCALE_9A ?= 3.125   # 9A one-third page, placed 2.5 in: 30 pt in R, 10.0 pt printed
FONT_MIN  ?= 10
FIG_SCALE ?= 1
# 'paper' reproduces the published panel's look (see --style in pca-stage.R).
STYLE     ?= current
LEGEND    ?= inside
YPAD      ?= 15

PLOT_SCRIPT := pca-stage.R
OUTDIR      := out/$(RUN)
STAMP       := $(OUTDIR)/.done

.PHONY: all paper download deps clean distclean print upload
.DEFAULT_GOAL := all

all: $(STAMP)

paper:
	$(MAKE) -B RUN=VZRJXQ-ruv-rbe STYLE=paper LEGEND=inside-genotype FONT_SCALE=$(PAPER_FONT_SCALE_1B) SHORT_STAGES=--short-stages
	$(MAKE) -B RUN=pca-hlh30-all-ruv STYLE=current LEGEND=inside-top FONT_SCALE=$(PAPER_FONT_SCALE_9A) SHORT_STAGES=--short-stages

download: $(COORDS) $(VAR)

$(COORDS):
	mkdir -p $(CACHE)
	aws s3 cp $(S3_COORDS) $@

$(VAR):
	mkdir -p $(CACHE)
	aws s3 cp $(S3_VAR) $@

$(STAMP): $(COORDS) $(VAR) $(PLOT_SCRIPT)
	mkdir -p $(OUTDIR)
	FONT_MIN=$(FONT_MIN) FIG_SCALE=$(FIG_SCALE) Rscript $(PLOT_SCRIPT) \
		--coords   $(COORDS) \
		--variance $(VAR) \
		--style    $(STYLE) \
		--legend   $(LEGEND) \
		--ypad     $(YPAD) \
		--font-scale $(FONT_SCALE) $(SHORT_STAGES) \
		--width    $(WIDTH) \
		--height   $(HEIGHT) \
		--outdir   $(OUTDIR)
	touch $@

deps:
	Rscript -e 'install.packages(c("ggplot2","ggforce","optparse","dplyr"), repos="https://cloud.r-project.org")'

S3_UPLOAD_DEST ?= $(S3_BUCKET)/posthoc_results/rz/pca-stage

# Syncs every run that has been built into out/<RUN>/ (no STAMP dep, so it
# uploads what's there regardless of which RUN was last built).
upload:
	aws s3 sync out $(S3_UPLOAD_DEST) --delete

clean:
	rm -rf out

distclean: clean
	rm -f $(COORDS) $(VAR)

print:
	@echo RUN=$(RUN)
	@echo S3_COORDS=$(S3_COORDS)
	@echo S3_VAR=$(S3_VAR)
	@echo OUTDIR=$(OUTDIR)
