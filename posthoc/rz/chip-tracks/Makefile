# HLH-30 ChIP panels for the hlh-30 paper (fig 8D/8E, fig S5A/S5B).
#
#   make            # render all four panels into out/
#   make tracks | schematic | motifs | lin
#
# Each script keeps its original canvas; FONT_SCALE_* raises the type so the
# smallest label prints at 10pt at the width the panel is placed in the
# author's layout (measured from Figures082826.pdf, see
# ../paper-figures/manifest.tsv):
#   FONT_SCALE = 10 / (smallest design size in pt * placed_width / canvas_width)
#   S5A metagene   canvas 5.4  placed 4.24  smallest 8.8pt  -> 1.447
#   S5B tracks     canvas 9.3  placed 7.09* smallest 5.5pt  -> 2.38   (*assumed)
#   8D  schematic  canvas 7.6  placed 5.80  smallest 8.5pt  -> 1.54
#   8E  motifs     canvas 4.1  placed 4.10  smallest 8pt    -> 1.25
#   S5B-opt lin    canvas 8.2  placed 7.09* smallest 6pt    -> 1.93   (*assumed)

OUTDIR ?= out
TRACK_GENES ?= ztf-28,tts-1,cpr-4
FONT_SCALE_METAGENE  ?= 1.447
FONT_SCALE_TRACKS    ?= 2.38
FONT_SCALE_SCHEMATIC ?= 1.54
FONT_SCALE_MOTIFS    ?= 1.25
FONT_SCALE_LIN       ?= 1.93

.PHONY: all tracks schematic motifs lin clean
.DEFAULT_GOAL := all

all: tracks schematic motifs lin

tracks:
	Rscript chip-tracks.R --genes "$(TRACK_GENES)" --outdir $(OUTDIR) \
	  --font-scale $(FONT_SCALE_METAGENE) --font-scale-tracks $(FONT_SCALE_TRACKS)

schematic:
	Rscript promoter-schematic.R --outdir $(OUTDIR) --font-scale $(FONT_SCALE_SCHEMATIC)

motifs:
	Rscript motif-counts.R --outdir $(OUTDIR) --font-scale $(FONT_SCALE_MOTIFS)

# Lin et al. 2018 Fig 5e layout: one averaged HLH-30 track per locus
lin:
	Rscript panel-e-lin.R --genes "$(TRACK_GENES)" --outdir $(OUTDIR) --font-scale $(FONT_SCALE_LIN)

clean:
	rm -f $(OUTDIR)/panelA_metagene.* $(OUTDIR)/panelE_tracks.* \
	      $(OUTDIR)/promoter_schematic.* $(OUTDIR)/motif_counts.* $(OUTDIR)/panelE_lin_style.*
