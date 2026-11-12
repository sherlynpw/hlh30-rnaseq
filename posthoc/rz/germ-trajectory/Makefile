# Germ-cell gene trajectories for fig S3, two gene sets side by side:
#   S3A  heatmaps, oocyte and germline genes (left) and sperm genes (right),
#        same row height and colour scale
#   S3B  oocyte and germline mean line, S3C sperm mean line, one y range
# Each panel is drawn at a canvas width and font scale that give 10pt in print:
#   S3A heatmaps     canvas 3.50 placed 3.50 gene names 8pt  -> 10/8             = 1.25
#   S3B/C mean lines canvas 4.4  placed 3.50 axis text 9.6pt -> 10/(9.6*0.795) = 1.31

FONT_SCALE_S3A ?= 1.25
FONT_SCALE_S3BC ?= 1.31
HM_S3A ?= --hm-width 3.5 --row-in 0.16 --zlim 2.1
YLIM_S3BC ?= -1.75,1.25
# the sperm set: the pseudopodium term holds the major sperm protein genes
SPERM_KW ?= pseudopod|\bsperm|spermat|male gamete
OOCYTE = --title "Oocyte and germline genes"
SPERM = --title "Sperm genes" --genes cache/sperm_genes.txt --prefix sperm

.PHONY: paper fetch genes stats clean
.DEFAULT_GOAL := paper

# corrected VST and sample sheet of hlh30-all-ruv, and the adult compare-DE GO
# tables (common filter, see ../compare-de) the gene set is built from
S3_BUCKET ?= s3://celegans-pipeline-rayandrew
fetch:
	mkdir -p cache
	aws s3 cp --only-show-errors $(S3_BUCKET)/results/hlh30-all-ruv/deseq2/vst_counts.csv cache/vst.csv
	aws s3 cp --only-show-errors $(S3_BUCKET)/results/hlh30-all-ruv/counts/sample_meta.tsv cache/meta.tsv
	aws s3 cp --only-show-errors --recursive $(S3_BUCKET)/results/VZRJXQ-ruv-rbe/enrichment/go/compare_de_tested_both/ \
	  cache/go/ --exclude "*" --include "*go_*.csv"

genes:
	Rscript germ-genes.R --go cache/go --out cache/germ_genes.txt --exclude '$(SPERM_KW)'
	Rscript germ-genes.R --go cache/go --out cache/sperm_genes.txt --kw '$(SPERM_KW)'

paper: genes
	Rscript germ-trajectory.R --font-scale $(FONT_SCALE_S3A) --outdir out-s3a $(HM_S3A) $(OOCYTE)
	Rscript germ-trajectory.R --font-scale $(FONT_SCALE_S3A) --outdir out-s3a $(HM_S3A) $(SPERM)
	Rscript germ-trajectory.R --font-scale $(FONT_SCALE_S3BC) --outdir out-s3b --width 4.4 --height 3.6 \
	  --ylim=$(YLIM_S3BC) $(OOCYTE)
	Rscript germ-trajectory.R --font-scale $(FONT_SCALE_S3BC) --outdir out-s3c --width 4.4 --height 3.6 \
	  --ylim=$(YLIM_S3BC) $(SPERM)

# WT vs hlh-30 at each stage, per gene set (camera and a per-sample score)
stats: genes
	Rscript germ-stats.R --out out-stats/germ_stats.tsv

clean:
	rm -rf out-s3a out-s3b out-s3c out-stats
