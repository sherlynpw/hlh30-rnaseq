# Fold-change quadrant plots, the dot-plot companion to the DEG Venns.
#   make paper      # four transitions into out/ (3.75 in) and out/third/ (2.5 in, fig 9I)
# Canvas 4.6 in: FONT_SCALE = 10 / (9.6 * placed / 4.6), rounded up to a whole point.
PAPER_FONT_SCALE ?= 1.3542      # 4A/4D/4G half page, placed 3.75 in: 13 pt in R, 10.6 pt printed
PAPER_FONT_SCALE_9I ?= 1.9792   # 9I one third, placed 2.5 in: 19 pt in R, 10.3 pt printed

.PHONY: paper clean
.DEFAULT_GOAL := paper

paper:
	Rscript fc-quadrant.R --font-scale $(PAPER_FONT_SCALE) --outdir out
	Rscript fc-quadrant.R --font-scale $(PAPER_FONT_SCALE_9I) --legend-ncol 1 --height 7.6 --outdir out/third

clean:
	rm -rf out
