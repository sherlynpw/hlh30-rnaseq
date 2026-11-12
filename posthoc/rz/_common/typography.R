# Shared typography floor for the posthoc/rz figure scripts. Sourced, not run
# directly.
#
# The journal requires every label in the assembled figure to print at >= 10pt.
# Panels are rendered here and then placed *scaled down* in the layout document,
# so a 10pt label in R prints smaller than 10pt on the page. FIG_SCALE carries
# that placement factor: a panel rendered 5.0in wide and placed 3.5in wide is
# FIG_SCALE = 0.70, and a label must be set at 10 / 0.70 = 14.3pt in R to land
# at 10pt in print.
#
#   FIG_SCALE=0.70 make -C volcano
#
# Per-panel factors: set FIG_SCALE in that subproject's Makefile.

FONT_MIN  <- as.numeric(Sys.getenv("FONT_MIN",  "10"))
FIG_SCALE <- as.numeric(Sys.getenv("FIG_SCALE", "1"))

stopifnot(is.finite(FONT_MIN), FONT_MIN > 0,
          is.finite(FIG_SCALE), FIG_SCALE > 0)

# Size in pt to request from R so the element prints at >= FONT_MIN once placed.
# Sizes below the floor are raised to it; sizes above keep their relative accent
# (a 12pt bold title stays visibly larger than 10pt body text).
fpt <- function(pt = FONT_MIN) pmax(pt, FONT_MIN) / FIG_SCALE

# geom_text()/geom_label()/geom_*_repel() take size in **mm**, not pt. Passing a
# pt value there silently renders ~2.8x too small -- always wrap it in fmm().
fmm <- function(pt = FONT_MIN) fpt(pt) / .pt

# Floor every text element a reader sees in the published panel. Provenance
# captions and subtitles are deliberately exempt: they are cropped or dropped
# during layout, and growing them reflows the panel for no printed benefit.
theme_floor <- function() {
  theme(
    axis.text    = element_text(size = fpt()),
    axis.text.x  = element_text(size = fpt()),
    axis.text.y  = element_text(size = fpt()),
    axis.title    = element_text(size = fpt()),
    legend.text  = element_text(size = fpt()),
    legend.title = element_text(size = fpt()),
    strip.text   = element_text(size = fpt())
  )
}
