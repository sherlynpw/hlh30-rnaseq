# Gene-list panels (Fig 1F, Fig S1A) from the L4+D3+D6 DESeq2 tables.
#   make     # out/fig-1f.*, out/fig-s1a.* and the gene tables
# Inputs come from the volcano cache; run that first.

.PHONY: all clean
all:
	Rscript gene-lists.R --outdir out

clean:
	rm -rf out
