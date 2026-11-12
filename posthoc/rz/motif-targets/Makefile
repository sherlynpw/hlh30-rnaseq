# What fraction of each DESeq2 DEG list are direct HLH-30 targets, defined by
# an HLH-30 E-box (~CAYRTG, JASPAR MA1449.1 PWM, p<1e-4) in the promoter
# -2000/+200 of the TSS, tested against the expressed-gene background.
#
#   make            # scan + stats + both figures for RUN=hlh30-all-ruv
#   make upload     # sync the dumbbell + summary to S3
#   make clean | distclean
#
# The DEG lists, DESeq2 background CSVs and GTF are reused from the sibling
# chip-tracks / chip-promoter caches, so no extra downloads are needed.

RUN       ?= hlh30-all-ruv
# fig S5C (canonical E-box only): 10 x 10 in canvas placed 5.0 in (two thirds of
# the page, from the draft's 4.79 in); smallest text 8.8pt x 2.2727 x 0.5 = 10pt in print
FONT_SCALE ?= 2.2727
WIDTH      ?= 10
HEIGHT     ?= 7.2
S3_BUCKET ?= s3://celegans-pipeline-rayandrew
GENE_DB   ?= ../../../gene_db.sqlite
CONTRASTS ?= tm1978_L1_vs_WT_L1,tm1978_L4_vs_WT_L4,tm1978_Day3_vs_WT_Day3,tm1978_Day6_vs_WT_Day6
WT_LABEL  ?= WT

DEDIR   := ../chip-promoter/cache/$(RUN)/de
DEGDIR  := ../chip-promoter/cache/$(RUN)/deg
COMPDIR := ../compare-de/out/$(RUN)/annotated
GTF     := cache/gtf.gtf
OUTDIR  := out/$(RUN)
STAMP   := $(OUTDIR)/.done

.PHONY: all clean distclean upload print
.DEFAULT_GOAL := all

all: $(STAMP)

$(STAMP): motif-targets.R render-targets.R direct-estimates.R $(GTF)
	Rscript motif-targets.R --run $(RUN) --gtf $(GTF) --gene-db $(GENE_DB) \
	  --de-dir $(DEDIR) --deg-dir $(DEGDIR) --compare-dir $(COMPDIR) \
	  --contrasts "$(CONTRASTS)"
	Rscript render-targets.R --summary $(OUTDIR)/motif_targets_summary.tsv \
	  --outdir $(OUTDIR) --wt-label $(WT_LABEL) --font-scale $(FONT_SCALE) \
	  --width $(WIDTH) --height $(HEIGHT)
	Rscript direct-estimates.R --summary $(OUTDIR)/motif_targets_summary.tsv \
	  --out $(OUTDIR)/direct_target_estimates.csv
	touch $@

# GTF is large; borrow it from the chip-tracks cache if present
$(GTF):
	mkdir -p cache
	cp ../chip-tracks/cache/gtf.gtf $(GTF)

S3_UPLOAD_DEST ?= $(S3_BUCKET)/posthoc_results/rz/motif-targets

# only the chosen figure (dumbbell) + the stats table
upload:
	aws s3 cp $(OUTDIR)/targets_dumbbell.png $(S3_UPLOAD_DEST)/targets_dumbbell.png
	aws s3 cp $(OUTDIR)/targets_dumbbell.pdf $(S3_UPLOAD_DEST)/targets_dumbbell.pdf
	aws s3 cp $(OUTDIR)/motif_targets_summary.tsv $(S3_UPLOAD_DEST)/motif_targets_summary.tsv

clean:
	rm -rf out

distclean: clean
	rm -rf cache

print:
	@echo RUN=$(RUN)
	@echo CONTRASTS=$(CONTRASTS)
	@echo WT_LABEL=$(WT_LABEL)
