# Percent of each timepoint's DEGs with experimental HLH-30 ChIP-seq binding
# near the TSS (ChIP-Atlas ce11), against the expressed-gene background.
#
#   make                 # DIST=1 and DIST=5 for RUN=hlh30-all-ruv
#   make DIST=1          # single window
#   make upload | clean | distclean

RUN       ?= hlh30-all-ruv
S3_BUCKET ?= s3://celegans-pipeline-rayandrew
CA_URL    ?= https://chip-atlas.dbcls.jp/data/ce11/target
DISTS     ?= 1 5
SCORE_MIN ?= 50
GENE_DB   ?= ../../../gene_db.sqlite

CONTRASTS ?= tm1978_L1_vs_WT_L1,tm1978_L4_vs_WT_L4,tm1978_Day3_vs_WT_Day3,tm1978_Day6_vs_WT_Day6

CACHE  := cache
DEGDIR := $(CACHE)/$(RUN)/deg
DEDIR  := $(CACHE)/$(RUN)/de
OUTDIR := out/$(RUN)
STAMP  := $(OUTDIR)/.done
PLOT_SCRIPT := chip-promoter.R

.PHONY: all deps clean distclean print upload
.DEFAULT_GOAL := all

all: $(STAMP)

$(STAMP): $(PLOT_SCRIPT)
	mkdir -p $(DEGDIR) $(DEDIR) $(OUTDIR)
	for d in $(DISTS); do \
	  [ -f $(CACHE)/hlh-30.$$d.tsv ] || curl -sSL -o $(CACHE)/hlh-30.$$d.tsv $(CA_URL)/hlh-30.$$d.tsv; \
	done
	for c in $$(echo $(CONTRASTS) | tr ',' ' '); do \
	  for k in de up down; do \
	    [ -f $(DEGDIR)/$${c}_genes_$$k.txt ] || \
	      aws s3 cp $(S3_BUCKET)/results/$(RUN)/filtered/filtered/$${c}_genes_$$k.txt $(DEGDIR)/ ; \
	  done; \
	  [ -f $(DEDIR)/$$c.csv ] || \
	    aws s3 cp $(S3_BUCKET)/results/$(RUN)/deseq2/results/$$c.csv $(DEDIR)/ ; \
	done
	for d in $(DISTS); do \
	  Rscript $(PLOT_SCRIPT) \
	    --chipatlas $(CACHE)/hlh-30.$$d.tsv \
	    --dist      $$d \
	    --deg-dir   $(DEGDIR) \
	    --de-dir    $(DEDIR) \
	    --contrasts "$(CONTRASTS)" \
	    --gene-db   $(GENE_DB) \
	    --score-min $(SCORE_MIN) \
	    --outdir    $(OUTDIR)/dist$$d; \
	done
	touch $@

deps:
	Rscript -e 'install.packages(c("optparse","dplyr","readr","ggplot2","RSQLite"), repos="https://cloud.r-project.org")'

S3_UPLOAD_DEST ?= $(S3_BUCKET)/posthoc_results/rz/chip-promoter

upload:
	aws s3 sync out $(S3_UPLOAD_DEST) --delete

clean:
	rm -rf out

distclean: clean
	rm -rf $(CACHE)

print:
	@echo RUN=$(RUN)
	@echo DISTS=$(DISTS)
	@echo SCORE_MIN=$(SCORE_MIN)
	@echo CONTRASTS=$(CONTRASTS)
