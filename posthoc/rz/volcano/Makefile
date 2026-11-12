# Volcano plots for mutant vs WT at each stage (L1, L4, Day3, Day6).
# Each per-stage DESeq2 contrast pools all temperatures
# (temperature is a covariate in the pipeline design, not a contrast term).
#
#   make                         # render all stages for default RUN
#   make PADJ_CUT=0.1 LFC_CUT=0.585
#   make RUN=<other-run>
#   make clean | distclean

S3_BUCKET ?= s3://celegans-pipeline-rayandrew
RUN       ?= hlh30-all-ruv
# Paper figure: WT vs hlh-30(tm1978) at L4, Day 3, Day 6 (see `paper`).
# L1 is still available with STAGES="L1 L4 Day3 Day6".
STAGES    ?= L4 Day3 Day6
MUT       ?= tm1978
WT        ?= WT
HIGHLIGHT ?= ztf-28,tts-1,cpr-4
PADJ_CUT  ?= 0.05
LFC_CUT   ?= 1
Y_MAX     ?= 30
LABEL_TOP ?= 0

# Paper panels keep the script's 5.0 x 4.6 canvas and raise FONT_SCALE so the
# axis text (0.8 x base 12) prints at 10pt at the placed width measured in
# Figures082826.pdf: 3.13 in for D0/D3/D6 (fig 1E), 2.50 in for L1 (fig 9D).
#   FONT_SCALE = 10 / (9.6 * placed / WIDTH)
# `make paper` renders both sets. Shrinking WIDTH to the placed size instead
# doubles every element relative to the plot.
WIDTH     ?= 5.0
HEIGHT    ?= 4.6
FONT_SCALE ?= 1
PAPER_FONT_SCALE    ?= 2.0833  # fig 1E one third of the page, placed 2.5 in, as 9D
PAPER_FONT_SCALE_L1 ?= 2.0833  # fig 9D one third of the page, placed 2.5 in (20 pt in R, 10.0 pt printed)
SHORT_STAGES ?=
SHORT_COUNTS ?=
STACK_COUNTS ?=
PLAIN_COUNTS ?=
LABEL_GAP    ?= 0.3
COUNT_SIZE   ?= 3.4
FONT_MIN  ?= 10
FIG_SCALE ?= 1

PLOT_SCRIPT := render-volcano.R
CACHE       := cache/$(RUN)
OUTDIR      := out/$(RUN)
STAMP       := $(OUTDIR)/.done

.PHONY: all paper deps clean distclean print upload
.DEFAULT_GOAL := all

all: $(STAMP)

# Fig 1E comes from the L4 + Day 3 + Day 6 model (VZRJXQ-ruv-rbe), the same
# run as the paper's DEG counts, Fig 1F, S1A and the clustering. Only the L1
# panel (Fig 9D) needs the four-stage model.
# No cap (Y_MAX=0): every dot at its real position, so no triangles; counts
# sit above the plot ("N down" flush left, "N up" flush right) and the stage
# name is inside, top left; labels sit further from their dots.
PAPER_OPTS = Y_MAX=0 PLAIN_COUNTS="--plain-counts --stage-inside" LABEL_GAP=1.2 SHORT_STAGES=--short-stages
paper:
	$(MAKE) -B RUN=VZRJXQ-ruv-rbe MUT=hlh-30 WT=N2 STAGES="L4 Day3 Day6" FONT_SCALE=$(PAPER_FONT_SCALE) $(PAPER_OPTS)
	$(MAKE) -B STAGES="L1" FONT_SCALE=$(PAPER_FONT_SCALE_L1) $(PAPER_OPTS)

$(STAMP): $(PLOT_SCRIPT)
	mkdir -p $(CACHE) $(OUTDIR)
	for s in $(STAGES); do \
	  f=$(MUT)_$${s}_vs_$(WT)_$${s}.csv; \
	  [ -f $(CACHE)/$$f ] || aws s3 cp $(S3_BUCKET)/results/$(RUN)/deseq2/results/$$f $(CACHE)/$$f; \
	done
	FONT_MIN=$(FONT_MIN) FIG_SCALE=$(FIG_SCALE) Rscript $(PLOT_SCRIPT) \
	  --cache    $(CACHE) \
	  --stages   "$(STAGES)" \
	  --mut      $(MUT) \
	  --wt       $(WT) \
	  --highlight "$(HIGHLIGHT)" \
	  --padj-cut $(PADJ_CUT) \
	  --lfc-cut  $(LFC_CUT) \
	  --y-max    $(Y_MAX) \
	  --label-top $(LABEL_TOP) \
	  --font-scale $(FONT_SCALE) $(SHORT_STAGES) $(SHORT_COUNTS) $(STACK_COUNTS) $(PLAIN_COUNTS) \
	  --label-gap $(LABEL_GAP) \
	  --count-size $(COUNT_SIZE) \
	  --width    $(WIDTH) \
	  --height   $(HEIGHT) \
	  --outdir   $(OUTDIR)
	touch $@

deps:
	Rscript -e 'install.packages(c("optparse","ggplot2","dplyr","readr","ggrepel","patchwork"), repos="https://cloud.r-project.org")'

S3_UPLOAD_DEST ?= $(S3_BUCKET)/posthoc_results/rz/volcano

upload:
	aws s3 sync out $(S3_UPLOAD_DEST) --delete

clean:
	rm -rf $(OUTDIR)

distclean: clean
	rm -rf $(CACHE)

print:
	@echo RUN=$(RUN)
	@echo STAGES=$(STAGES)
	@echo HIGHLIGHT=$(HIGHLIGHT)
	@echo PADJ_CUT=$(PADJ_CUT)
	@echo LFC_CUT=$(LFC_CUT)
	@echo Y_MAX=$(Y_MAX)
	@echo LABEL_TOP=$(LABEL_TOP)
	@echo OUTDIR=$(OUTDIR)
