# Collects the rendered panels from the sibling subprojects under their figure
# names (fig-1b.pdf, fig-1c.pdf, ...) and syncs them to S3, so a collaborator
# can find a panel by its label in the manuscript rather than by which script
# happened to produce it.
#
#   make            # stage out/fig-*.pdf|png from manifest.tsv
#   make new        # stage the new analyses (S7 to S11) into out/
#   make go-full    # stage out/fig-*-full.xlsx (every significant GO term)
#   make upload     # stage both, then sync out/ to S3 (PDFs, tables and docx; PNGs stay local)
#   make check      # list manifest rows whose source is missing
#   make measure    # printed point size of every staged panel (needs uv)
#   make sheets     # out/ALL-figures.pdf, for local review (not uploaded)
#
# Staged PDFs have every font embedded (see EMBED).
#
# manifest.tsv is the single source of truth: panel <TAB> source <TAB> size.
# Add a row when a new panel is rendered. Sources are PDFs; the matching PNG is
# copied too when it exists.

S3_BUCKET ?= s3://celegans-pipeline-rayandrew
S3_DEST   ?= $(S3_BUCKET)/posthoc_results/rz/paper-figures

MANIFEST := manifest.tsv
# R's pdf() names Helvetica without embedding it, so a reader's PDF program
# swaps in its own font (often Arial). Ghostscript rewrites each staged PDF
# with every font embedded (Helvetica as its metric clone Nimbus Sans; the
# empty NeverEmbed list is what makes gs embed the base-14 fonts too); the
# drawing and any raster images are kept as they are.
EMBED = gs -q -dNOPAUSE -dBATCH -dSAFER -sDEVICE=pdfwrite -dEmbedAllFonts=true \
	-dSubsetFonts=true -dAutoRotatePages=/None -dDownsampleColorImages=false \
	-dDownsampleGrayImages=false -dDownsampleMonoImages=false \
	-dAutoFilterColorImages=false -dColorImageFilter=/FlateEncode
OUTDIR   := out

NEW      := new-analyses.tsv

.PHONY: all new go-full upload check measure sheets clean print
.DEFAULT_GOAL := all

all:
	@mkdir -p $(OUTDIR)
	@# clears only the manifest's own panels, so the new analyses staged in the
	@# same folder by `make new` are kept
	@n=0; miss=0; TAB=$$(printf '\t'); \
	while IFS="$$TAB" read -r panel src size note; do \
	  case "$$panel" in ''|\#*) continue;; esac; \
	  rm -f "$(OUTDIR)/$$panel.pdf" "$(OUTDIR)/$$panel.png"; \
	  if [ -f "../$$src" ]; then \
	    $(EMBED) -sOutputFile="$(OUTDIR)/$$panel.pdf" -c "<< /NeverEmbed [ ] >> setdistillerparams" -f "../$$src" || exit 1; \
	    png="../$${src%.pdf}.png"; \
	    [ -f "$$png" ] && cp "$$png" "$(OUTDIR)/$$panel.png"; \
	    n=$$((n+1)); \
	  else \
	    echo "  MISSING $$panel <- $$src"; miss=$$((miss+1)); \
	  fi; \
	done < $(MANIFEST); \
	echo "staged $$n panel(s) into $(OUTDIR)/; $$miss missing"

check:
	@TAB=$$(printf '\t'); \
	while IFS="$$TAB" read -r panel src size note; do \
	  case "$$panel" in ''|\#*) continue;; esac; \
	  [ -f "../$$src" ] || echo "MISSING $$panel <- $$src"; \
	done < $(MANIFEST); echo "check done"

# Paper versions (10pt at their placed width, no notes) go in out/ next to the
# panels.
# PDFs get the same font embedding as the panels.
new:
	@$(call stage_list,$(NEW),$(OUTDIR),)

# $(1) = list (name <TAB> source <TAB> note [<TAB> pages, e.g. 5-8]), $(2) = folder,
# $(3) = "own" when the folder holds only this list: then its files are cleared
# first and README.txt is written from the notes. Otherwise only the listed
# names are replaced.
define stage_list
mkdir -p $(2); [ "$(3)" = own ] && find $(2) -maxdepth 1 -type f -delete; \
n=0; miss=0; TAB=$$(printf '\t'); \
while IFS="$$TAB" read -r name src note pages; do \
  case "$$name" in ''|\#*) continue;; esac; \
  rm -f "$(2)/$$name".*; \
  if [ -f "../$$src" ]; then \
    ext="$${src##*.}"; \
    if [ "$$ext" = pdf ]; then \
      rng=""; [ -n "$$pages" ] && rng="-dFirstPage=$${pages%-*} -dLastPage=$${pages#*-}"; \
      $(EMBED) $$rng -sOutputFile="$(2)/$$name.pdf" -c "<< /NeverEmbed [ ] >> setdistillerparams" -f "../$$src" || exit 1; \
    else cp "../$$src" "$(2)/$$name.$$ext"; fi; \
    png="../$${src%.*}.png"; \
    [ "$$ext" = pdf ] && [ -f "$$png" ] && cp "$$png" "$(2)/$$name.png"; \
    n=$$((n+1)); \
  else \
    echo "  MISSING $$name <- $$src"; miss=$$((miss+1)); \
  fi; \
done < $(1); \
[ "$(3)" = own ] && grep -v '^#' $(1) | awk -F'\t' 'NF>=3 {printf "- %s: %s\n", $$1, $$3}' > $(2)/README.txt; \
echo "staged $$n file(s) into $(2)/; $$miss missing"
endef

# One sheet per GO panel with every significant term, the panel's terms bold.
go-full:
	$(MAKE) -C ../go-plot full
	@mkdir -p $(OUTDIR)
	@rm -f $(OUTDIR)/fig-*-full.xlsx
	cp ../go-plot/out-full/fig-*-full.xlsx $(OUTDIR)/

# Reads each panel's text matrices and scales by placed / canvas width.
measure:
	uv run --with pikepdf python measure-pt.py

# One review PDF of the manifest panels, for local use.
sheets: all new
	Rscript contact-sheet.R --indir $(OUTDIR) --out $(OUTDIR)/ALL-figures.pdf

# No --delete: this prefix is shared output, never prune what we did not stage.
upload: all new go-full
	aws s3 sync $(OUTDIR) $(S3_DEST) --exclude '*ALL-figures.pdf' --exclude '*README.txt' --exclude '*.png'

clean:
	rm -rf $(OUTDIR)

print:
	@echo "dest: $(S3_DEST)"; echo "panels:"; grep -vc '^#' $(MANIFEST)
