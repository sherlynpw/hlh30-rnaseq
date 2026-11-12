# Re-render the genotype-PCA from a pipeline run's exploration/pca outputs.
# Same family as posthoc/rz/pca-stage but groups ellipses by genotype x stage
# (4 ellipses for the no_l4 run: N2 D3, N2 D6, hlh-30 D3, hlh-30 D6).
#
#   make                          # download + plot for RUN=VZRJXQ_no_l4-ruv-rbe
#   make RUN=<other-run>
#   make upload                   # sync out/ back to S3
#   make clean | distclean

RUN       ?= VZRJXQ_no_l4-ruv-rbe
S3_BUCKET ?= s3://celegans-pipeline-rayandrew

S3_PCA    := $(S3_BUCKET)/results/$(RUN)/exploration/pca
S3_COORDS := $(S3_PCA)/pca_coordinates.csv
S3_VAR    := $(S3_PCA)/variance_explained.txt

# Namespaced by RUN: a bare filename meant switching RUN silently reused the
# previously downloaded CSV, so out/<RUN>/ could hold another run's plot.
CACHE  := cache/$(RUN)
COORDS := $(CACHE)/pca_coordinates.csv
VAR    := $(CACHE)/variance_explained.txt

# Fig 1D is placed at 3.75 x 2.75 in (95.2 x 69.8 mm), measured with pdfimages.
# Keep the 7.5 in canvas; `make paper` raises the type so axis text
# (0.8 x base 12) prints at 10pt: FONT_SCALE = 10 / (9.6 * 3.75 / 7.5).
WIDTH     ?= 7.5
HEIGHT    ?= 5.5
FONT_SCALE ?= 1
GROUP_SHAPE ?= hull
PAPER_FONT_SCALE ?= 2.0833  # 1D half page, placed 3.75 in (20 pt in R, 10.0 pt printed)
SHORT_STAGES ?=
FONT_MIN  ?= 10
FIG_SCALE ?= 1
LEGEND    ?= inside

PLOT_SCRIPT := pca-genotype.R
OUTDIR      := out/$(RUN)
STAMP       := $(OUTDIR)/.done

.PHONY: all paper download deps clean distclean print upload
.DEFAULT_GOAL := all

all: $(STAMP)

paper:
	$(MAKE) -B FONT_SCALE=$(PAPER_FONT_SCALE) SHORT_STAGES=--short-stages LEGEND=inside-h GROUP_SHAPE=mark
	$(MAKE) -B FONT_SCALE=$(PAPER_FONT_SCALE) SHORT_STAGES=--short-stages LEGEND=inside-h GROUP_SHAPE=mark \
		COLOR_BY=cohort OUTNAME=pca_cohort

# fig S1C: the fig 1D PCA coloured by adult cohort (20 C or room temperature)
COLOR_BY ?= stage
OUTNAME  ?= pca_genotype

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
		--legend   $(LEGEND) \
		--group-shape $(GROUP_SHAPE) \
		--font-scale $(FONT_SCALE) $(SHORT_STAGES) \
		--width    $(WIDTH) \
		--height   $(HEIGHT) \
		--color-by $(COLOR_BY) --outname $(OUTNAME) \
		--outdir   $(OUTDIR)
	touch $@

deps:
	Rscript -e 'install.packages(c("ggplot2","ggforce","optparse","dplyr"), repos="https://cloud.r-project.org")'

S3_UPLOAD_DEST ?= $(S3_BUCKET)/posthoc_results/rz/pca-genotype

upload:
	aws s3 sync out $(S3_UPLOAD_DEST) --delete

clean:
	rm -rf out

distclean: clean
	rm -f $(COORDS) $(VAR)

print:
	@echo RUN=$(RUN)
	@echo S3_COORDS=$(S3_COORDS)
	@echo OUTDIR=$(OUTDIR)
