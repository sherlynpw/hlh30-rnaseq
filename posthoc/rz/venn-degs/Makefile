# DEG Venn diagrams, WT vs hlh-30(tm1978) per transition (see README.md).
#
#   make paper      # figs 4A/4D/4G into out/, fig 9I into out9/
#
# All four render on the script's 4.6 x 3.9 in canvas; FONT_SCALE raises the
# type so the smallest label (set titles, 2.7 mm = 7.7 pt) prints at 10pt at
# the placed width: 3.6 in for all four (9I was placed 2.48 in, but its
# near-concentric circles cannot hold 10pt side labels there).
#   FONT_SCALE = 10 / (7.7 * placed / 4.6)

PAPER_FONT_SCALE    ?= 1.66
PAPER_FONT_SCALE_9I ?= 1.66

.PHONY: paper clean
.DEFAULT_GOAL := paper

paper:
	Rscript venn-degs.R --short-stages --font-scale $(PAPER_FONT_SCALE) --outdir out
	Rscript venn-degs.R --short-stages --font-scale $(PAPER_FONT_SCALE_9I) --outdir out9

clean:
	rm -rf out out9
