#!/usr/bin/env Rscript
# Build a single review PDF from the staged panels: one page per panel, each
# labelled with its figure name, so the whole set can be checked in one file
# instead of opening 30-odd PDFs. Reads the PNGs (not the PDFs) because they
# are already the rendered result and embed without a PDF toolchain.
#
#   Rscript contact-sheet.R --indir out --out out/ALL-figures.pdf

suppressPackageStartupMessages({
  library(optparse); library(png); library(grid)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--indir", type = "character", default = "out"),
  make_option("--out",   type = "character", default = "out/ALL-figures.pdf"),
  make_option("--manifest", type = "character", default = "manifest.tsv")
)))

# Panels are drawn at their TRUE placed size (column 3 = placed width in
# inches), not scaled to fill the page. Fitting each panel to the page makes a
# wide panel's 10pt type look smaller than a tall panel's identical 10pt, which
# defeats the point of a sheet you review type size on.
read_manifest <- function(path) {
  ln <- readLines(path)
  # the header line itself starts with "#", so the comment filter already
  # removed it -- dropping another line here would lose the first panel
  ln <- ln[!grepl("^\\s*#", ln) & nzchar(trimws(ln))]
  f  <- strsplit(ln, "\t")
  data.frame(panel = vapply(f, `[`, character(1), 1),
             placed_w = suppressWarnings(as.numeric(vapply(f, `[`, character(1), 3))),
             stringsAsFactors = FALSE)
}

mf     <- read_manifest(opt$manifest)
panels <- mf$panel
files  <- file.path(opt$indir, paste0(panels, ".png"))
keep   <- file.exists(files)
if (any(!keep)) message("[skip] no PNG for: ", paste(panels[!keep], collapse = ", "))
panels <- panels[keep]; files <- files[keep]; placed_w <- mf$placed_w[keep]
stopifnot(length(files) > 0)

# A4 landscape gives wide panels room; each page is one panel, scaled to fit
# while preserving its aspect ratio, with the figure name across the top.
pdf(opt$out, width = 11.69, height = 8.27, onefile = TRUE)
for (i in seq_along(files)) {
  img <- readPNG(files[i])
  ar  <- dim(img)[1] / dim(img)[2]               # height / width
  grid.newpage()
  grid.text(panels[i], x = 0.02, y = 0.965, just = c("left", "top"),
            gp = gpar(fontsize = 20, fontface = "bold"))
  grid.text(sprintf("%d / %d", i, length(files)), x = 0.98, y = 0.965,
            just = c("right", "top"), gp = gpar(fontsize = 12, col = "grey45"))
  # true size: placed width in inches, as a fraction of the 11.69in page
  pw <- placed_w[i]
  shrunk <- FALSE
  w_in <- pw; h_in <- pw * ar
  if (h_in > 7.0) { h_in <- 7.0; w_in <- h_in / ar; shrunk <- TRUE }  # tall panels
  grid.raster(img, x = 0.5, y = 0.45,
              width = unit(w_in, "in"), height = unit(h_in, "in"))
  grid.text(sprintf("placed %.2f in (%.0f mm) wide%s", pw, pw * 25.4,
                    if (shrunk) "  - shown reduced to fit page height" else "  - shown at true size"),
            x = 0.02, y = 0.925, just = c("left", "top"),
            gp = gpar(fontsize = 9, col = "grey45"))
}
invisible(dev.off())
message("wrote ", opt$out, "  (", length(files), " panels)")
