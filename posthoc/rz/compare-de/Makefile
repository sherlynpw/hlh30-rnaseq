# Compare-DE (WT vs hlh-30 per stage transition), with GO on each gene list.
#
# compare_de.R on the run's DESeq2 tables, annotated up/down tables, then GO
# on each gene list in the pipeline's R container with the run's own universe.
# Only genes with a padj in BOTH contrasts are compared (--tested_in_both), so
# a gene is not "exclusive" to one genotype just because DESeq2 filtered it
# out of the other genotype's test. The Venns and quadrant plots use the same
# rule. `make upload` puts the lists under posthoc_results/rz/compare-de/ and
# the GO tables under results/<run>/enrichment/go/$(GO_DIR)/, which go-plot
# reads.
#
#   make RUN=hlh30-all-ruv     # L1 -> L4 (the pipeline run has compare_de off)
#   make RUN=VZRJXQ-ruv-rbe    # L4 -> D3, L4 -> D6, D3 -> D6
#   make upload RUN=...

RUN  ?= hlh30-all-ruv
B    ?= s3://celegans-pipeline-rayandrew
IMG  ?= rayandrew/celegans-r-analysis:0.11.0
ifeq ($(RUN),VZRJXQ-ruv-rbe)
PAIR := hlh-30_Day3_vs_hlh-30_L4__X__N2_Day3_vs_N2_L4,hlh-30_Day6_vs_hlh-30_L4__X__N2_Day6_vs_N2_L4,hlh-30_Day6_vs_hlh-30_Day3__X__N2_Day6_vs_N2_Day3
AGES := L4,Day3,Day6
GENO := N2,hlh-30
# the pipeline's own compare_de GO (no common filter) stays where it is
GO_DIR := compare_de_tested_both
else
PAIR := WT_L4_vs_WT_L1__X__tm1978_L4_vs_tm1978_L1
AGES := L1,L4,Day3,Day6
GENO := WT,tm1978
GO_DIR := compare_de
endif
REPO := $(abspath ../../..)
C    := cache/$(RUN)
O    := out/$(RUN)

.PHONY: all de compare annotate go upload clean
all: go

de: $(C)/de/.done
$(C)/de/.done:
	mkdir -p $(C)/de
	aws s3 cp $(B)/results/$(RUN)/deseq2/results/ $(C)/de/ --recursive --exclude "*" --include "*.csv" --only-show-errors
	aws s3 cp $(B)/results/$(RUN)/enrichment/background.txt $(C)/background.txt --only-show-errors
	touch $@

compare: de
	mkdir -p $(O)
	Rscript $(REPO)/bin/compare_de.R --input_dir $(C)/de --pairs "$(PAIR)" \
	  --age_order "$(AGES)" --geno_order "$(GENO)" --tested_in_both TRUE \
	  --padj_threshold 0.05 --lfc_threshold 1 --outdir $(O)
	for f in $(O)/*_compared_genes.csv; do Rscript split-differential.R --compared $$f || exit 1; done

annotate: compare
	for f in $(O)/*_compared_genes.csv; do \
	  Rscript annotate-lists.R --compared $$f --gene-db $(REPO)/gene_db.sqlite --outdir $(O)/annotated || exit 1; \
	done

# GO in the pipeline container, as the pipeline's GO step runs
go: annotate
	rm -rf $(C)/go-input $(O)/go && mkdir -p $(C)/go-input $(O)/go
	cp $(O)/*_genes_*.txt $(C)/go-input/
	cp $(REPO)/gene_db.sqlite $(C)/gene_db.sqlite
	docker run --rm --user $$(id -u):$$(id -g) -w /repo -v $(REPO):/repo:ro -v $(abspath $(C)):/work \
	  -v $(abspath $(O)/go):/out $(IMG) Rscript bin/go_enrichment.R --input /work/go-input \
	  --gene_db /work/gene_db.sqlite --universe /work/background.txt \
	  --pvalue 0.05 --qvalue 0.2 --cpus 6 --outdir /out

upload:
	aws s3 cp $(O)/annotated/ $(B)/posthoc_results/rz/compare-de/$(RUN)/annotated/ --recursive --exclude "*" --include "*.tsv" --only-show-errors
	aws s3 cp $(O)/ $(B)/posthoc_results/rz/compare-de/$(RUN)/gene_lists/ --recursive --exclude "*" --include "*_genes_*.txt" --exclude "go/*" --exclude "annotated/*" --only-show-errors
	aws s3 cp $(O)/compare_de_report.csv $(B)/posthoc_results/rz/compare-de/$(RUN)/compare_de_report.csv --only-show-errors
	aws s3 cp $(O)/go/ $(B)/results/$(RUN)/enrichment/go/$(GO_DIR)/ --recursive --exclude "*" --include "*.csv" --include "*.pdf" --only-show-errors

clean:
	rm -rf cache out
