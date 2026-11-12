#!/usr/bin/env Rscript
# Area-proportional 2-set Venn diagrams of WT vs hlh-30(tm1978) DEG sets, one per
# developmental transition (L1->L4, L4->Day3, L4->Day6, Day3->Day6). Exclusive
# regions are split up/down; the overlap is split up/down too.
#
# Geometry: for two circles there is an exact area-proportional layout - circle
# areas proportional to set sizes, centre distance solved so the lens (overlap)
# area matches the shared count. No eulerr dependency (2 sets solve analytically).
#
# Two overlap definitions are rendered (see README):
#   shared        overlap = concordant "shared" only (DE in both, |dlog2FC|<1),
#                 matching the compare-DE `shared` category. Differential genes
#                 (DE in both, |dlog2FC|>=1) are noted, not drawn in the overlap.
#   intersection  overlap = ALL genes DE in both genotypes (shared+differential),
#                 the strict set-theoretic intersection.
#
# All counts are computed here from the DESeq2 result CSVs so the figure and the
# README's numbers come from one place. DEG = padj<0.05 & |log2FC|>=1.

suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(ggforce); library(dplyr); library(readr); library(patchwork)
})

# >= 10pt in print: see _common/typography.R for FIG_SCALE / FONT_MIN.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "typography.R"))

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--bucket", type = "character", default = "s3://celegans-pipeline-rayandrew"),
  make_option("--cache",  type = "character", default = "cache"),
  make_option("--outdir", type = "character", default = "out"),
  make_option("--padj",   type = "double",    default = 0.05),
  make_option("--lfc",    type = "double",    default = 1.0),
  make_option("--short-stages", action = "store_true", default = FALSE, dest = "short_stages",
              help = "Use the short stage labels L1/L4/D3/D6, matching the rest of the figure set."),
  make_option("--l4-name", type = "character", default = "L4", dest = "l4_name",
              help = "Short label for the L4 stage with --short-stages (paper panels: L4)."),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size. Render at the ORIGINAL canvas and raise this."),
  make_option("--width",  type = "double",    default = 4.6,
              help = "Panel width in inches. Set to the placed width for 1:1 layout."),
  make_option("--height", type = "double",    default = 3.9)
)))
dir.create(opt$cache, showWarnings = FALSE, recursive = TRUE)

# standard genotype palette (WT blue / hlh-30 orange), ztf-28 fill-only style
WTc <- PAL$wt; MUTc <- PAL$hlh30
UPc <- PAL$up; DNc <- PAL$down; SHc <- "grey15"; DIFc <- "grey35"

# --- transitions: WT contrast + hlh-30 contrast, with orientation (flip) so that
#     "up" always means higher at the LATER stage of the transition -------------
TR <- list(
  list(id="A", title=if (opt$short_stages) paste("L1 →", opt$l4_name) else "L1 → L4", run="hlh30-all-ruv",
       wt="WT_L1_vs_WT_L4", wt_flip=TRUE, mut="tm1978_L4_vs_tm1978_L1", mut_flip=FALSE),
  list(id="B", title=if (opt$short_stages) paste(opt$l4_name, "→ D3") else "L4 → Day 3", run="VZRJXQ-ruv-rbe",
       wt="N2_Day3_vs_N2_L4", wt_flip=FALSE, mut="hlh-30_L4_vs_hlh-30_Day3", mut_flip=TRUE),
  list(id="C", title=if (opt$short_stages) paste(opt$l4_name, "→ D6") else "L4 → Day 6", run="VZRJXQ-ruv-rbe",
       wt="N2_Day6_vs_N2_L4", wt_flip=FALSE, mut="hlh-30_L4_vs_hlh-30_Day6", mut_flip=TRUE),
  list(id="D", title=if (opt$short_stages) "D3 → D6" else "Day 3 → Day 6", run="VZRJXQ-ruv-rbe",
       wt="N2_Day6_vs_N2_Day3", wt_flip=FALSE, mut="hlh-30_Day6_vs_hlh-30_Day3", mut_flip=FALSE)
)

fetch <- function(run, name) {
  f <- file.path(opt$cache, paste0(run, "__", name, ".csv"))
  if (!file.exists(f)) {
    s3 <- sprintf("%s/results/%s/deseq2/results/%s.csv", opt$bucket, run, name)
    message("[download] ", s3)
    if (system2("aws", c("s3", "cp", s3, f)) != 0) stop("download failed: ", s3)
  }
  f
}
tested <- function(run, name) {
  d <- read.csv(fetch(run, name), stringsAsFactors = FALSE)
  d$gene_id[!is.na(d$padj) & !is.na(d$log2FoldChange)]
}
# DEGs among `keep` only. Both genotypes are compared on the genes that DESeq2
# tested in both (a padj in both contrasts), so a gene is not "WT only" just
# because independent filtering dropped it from the hlh-30 test.
deg <- function(run, name, flip, keep) {
  d <- read.csv(fetch(run, name), stringsAsFactors = FALSE)
  d <- d[!is.na(d$padj) & !is.na(d$log2FoldChange) & d$gene_id %in% keep, ]
  if (flip) d$log2FoldChange <- -d$log2FoldChange
  d <- d[d$padj < opt$padj & abs(d$log2FoldChange) >= opt$lfc, c("gene_id", "log2FoldChange")]
  setNames(d$log2FoldChange, d$gene_id)
}

# counts for one transition; returns a one-row tibble with everything the figure
# and the README need.
counts_for <- function(t) {
  keep <- intersect(tested(t$run, t$wt), tested(t$run, t$mut))
  wt <- deg(t$run, t$wt, t$wt_flip, keep); mut <- deg(t$run, t$mut, t$mut_flip, keep)
  wg <- names(wt); mg <- names(mut)
  wo <- setdiff(wg, mg); mo <- setdiff(mg, wg); sh <- intersect(wg, mg)
  dl <- wt[sh] - mut[sh]
  shared_ids <- sh[abs(dl) < opt$lfc]; diff_ids <- sh[abs(dl) >= opt$lfc]
  tibble(
    id=t$id, title=t$title, run=t$run, wt_contrast=t$wt, mut_contrast=t$mut,
    wt_total=length(wg), mut_total=length(mg),
    wt_only_down=sum(wt[wo]<0), wt_only_up=sum(wt[wo]>0),
    mut_only_down=sum(mut[mo]<0), mut_only_up=sum(mut[mo]>0),
    n_intersect=length(sh),
    shared=length(shared_ids), shared_up=sum(wt[shared_ids]>0), shared_down=sum(wt[shared_ids]<0),
    differential=length(diff_ids),
    inter_up=sum(wt[sh]>0 & mut[sh]>0), inter_down=sum(wt[sh]<0 & mut[sh]<0),
    inter_mixed=sum(sign(wt[sh])!=sign(mut[sh]))
  )
}

# ------------------------------------------------------------ geometry ---------
# lens (intersection) area of two circles radius r1,r2, centre distance d
lens_area <- function(d, r1, r2) {
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  a1 <- r1^2 * acos((d^2 + r1^2 - r2^2) / (2*d*r1))
  a2 <- r2^2 * acos((d^2 + r2^2 - r1^2) / (2*d*r2))
  a3 <- 0.5 * sqrt((-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2))
  a1 + a2 - a3
}
solve_d <- function(r1, r2, target) {                 # centre distance for a given lens area
  target <- min(target, pi*min(r1,r2)^2 * 0.999)      # guard: can't exceed smaller circle
  if (target <= 0) return(r1 + r2)
  uniroot(function(d) lens_area(d, r1, r2) - target,
          lower = abs(r1-r2) + 1e-6, upper = r1 + r2 - 1e-6)$root
}

txt <- function(x, y, label, colour, size, face="plain", ...)
  annotate("text", x=x, y=y, label=label, colour=colour, size=size, fontface=face, ...)

# one Venn. mode: "shared" or "intersection".
draw_venn <- function(r, mode) {
  ov <- if (mode == "shared") r$shared else r$n_intersect
  ov_up   <- if (mode == "shared") r$shared_up   else r$inter_up
  ov_down <- if (mode == "shared") r$shared_down else r$inter_down
  wt_tot  <- r$wt_only_down + r$wt_only_up + ov
  mut_tot <- r$mut_only_down + r$mut_only_up + ov

  r1 <- sqrt(wt_tot/pi); r2 <- sqrt(mut_tot/pi)
  d  <- solve_d(r1, r2, ov)
  cx1 <- -d/2; cx2 <- d/2                              # WT left, hlh-30 right
  circ <- data.frame(x0=c(cx1,cx2), y0=0, r=c(r1,r2), set=c("WT","MUT"))

  # label anchors: centre of each exclusive crescent, and centre of the lens.
  # Using the true crescent/lens x-spans keeps labels out of the overlap even
  # when the overlap is huge (near-concentric circles).
  xw <- ((cx1 - r1) + (cx2 - r2)) / 2                  # WT-only crescent centre
  xm <- ((cx1 + r1) + (cx2 + r2)) / 2                  # hlh-30-only crescent centre
  xo <- ((cx2 - r2) + (cx1 + r1)) / 2                  # lens (overlap) centre
  # Every offset is a fraction of the largest radius, so with coord_fixed the
  # gaps are the same fraction of the panel whatever the gene counts; the title
  # gap is in pt, so it is identical across panels.
  sz <- max(r1, r2)                                    # text offset scale
  o  <- sz*0.16
  # A thin exclusive crescent (near-concentric circles, e.g. L1 -> L4) cannot
  # hold "1130 Down": put that side's counts just outside its circle instead.
  w_wt  <- (cx2 - r2) - (cx1 - r1)
  w_mut <- (cx2 + r2) - (cx1 + r1)
  thin  <- 0.25 * sz
  wt_out  <- w_wt  < thin
  mut_out <- w_mut < thin
  any_out <- wt_out || mut_out
  # Near-concentric circles (L1 -> L4, fig 9I): no crescent can hold a count,
  # so each set's exclusive counts sit BESIDE its circle, WT on the left
  # (right-aligned to the circle) and hlh-30 on the right. The x-range is
  # padded by the width of the longest such label so nothing is clipped.
  num <- 3.3 * opt$font_scale
  hw   <- if (any_out) 1 else 0.5
  hm   <- if (any_out) 0 else 0.5
  side_chars <- max(nchar(c(paste(r$wt_only_down, "Down"), paste(r$wt_only_up, "Up"),
                            paste(r$mut_only_down, "Down"), paste(r$mut_only_up, "Up"))))
  # One text line in data units: the panel is width-limited under coord_fixed,
  # so data-per-pt follows the x-range over the drawable width. Row offsets in
  # lines keep the stack legible at any font scale or gene count.
  draw_w <- opt$width * 72 - 16
  draw_h <- opt$height * 72 - 12 - 14 * opt$font_scale   # minus the title
  pad <- sz*0.45
  for (i in 1:4) {                                       # pad and lh depend on each other
    xr  <- (cx2 + r2 + pad) - (cx1 - r1 - pad)
    lh  <- 1.15 * num * .pt * xr / draw_w
    # ~0.55 em per character: room for the side labels plus a small gap
    if (any_out) pad <- max(pad, (side_chars * 0.55 + 0.6) * lh / 1.15)
    ylim <- c(-sz - 0.2*lh, sz + 1.1*lh)
    # keep the panel width-limited, or the line unit above is wrong
    need <- diff(ylim) * draw_w / draw_h
    if (xr < need) pad <- pad + (need - xr) / 2
  }
  y_dn <-  0.5*lh
  y_up <- -0.5*lh
  ytop <- sz + 0.35*lh
  # crescent labels sit a little outward of the crescent centre so the wider
  # "N Down" in the lens never meets them
  xw_n <- if (any_out) cx1 - r1 - 0.3*lh else xw - 0.08*sz
  xm_n <- if (any_out) cx2 + r2 + 0.3*lh else xm + 0.08*sz
  # set titles: over the crescents, or over the circle edges when there are none
  xw_t <- if (any_out) cx1 - r1 - pad else xw
  xm_t <- if (any_out) cx2 + r2 + pad else xm
  hw_t <- if (any_out) 0 else 0.5
  hm_t <- if (any_out) 1 else 0.5

  g <- ggplot() +
    geom_circle(data=circ, aes(x0=x0, y0=y0, r=r, fill=set), colour=NA, alpha=0.28) +
    scale_fill_manual(values=c(WT=WTc, MUT=MUTc), guide="none") +
    # set titles over each set's own region, coloured to match, same height
    txt(xw_t, ytop, "WT DEGs", WTc, 2.7 * opt$font_scale, "bold", vjust=0, hjust=hw_t) +
    # gene name italic, "DEGs" roman
    annotate("text", x = xm_t, y = ytop, label = 'bold(italic("hlh-30(tm1978)")~DEGs)', parse = TRUE,
             colour = MUTc, size = 2.7 * opt$font_scale, vjust = 0, hjust = hm_t) +
    # WT-only  (down over up)
    txt(xw_n, y_dn, paste(r$wt_only_down, "Down"), DNc, num, hjust=hw) +
    txt(xw_n, y_up, paste(r$wt_only_up,   "Up"),   UPc, num, hjust=hw) +
    # hlh-30-only
    txt(xm_n, y_dn, paste(r$mut_only_down, "Down"), DNc, num, hjust=hm) +
    txt(xm_n, y_up, paste(r$mut_only_up,   "Up"),   UPc, num, hjust=hm) +
    # overlap (compact up/down split)
    txt(xo,  1.5*lh, if (mode=="shared") "Shared" else "DE in both", SHc, 3.1 * opt$font_scale, "bold") +
    txt(xo,  0.5*lh, paste(ov_down, "Down"), DNc, num) +
    txt(xo, -0.5*lh, paste(ov_up,   "Up"),   UPc, num) +
    labs(title = r$title) +
    coord_fixed(clip="off", xlim=c(cx1 - r1 - pad, cx2 + r2 + pad), ylim=ylim, expand=FALSE) +
    theme_void(base_size=12 * opt$font_scale) +
    theme(plot.title = element_text(face="bold", size=10 * opt$font_scale, hjust=0.5,
                                    margin=margin(b=4 * opt$font_scale, unit="pt")),
          plot.margin = margin(6, 8, 6, 8))
  # in shared mode, note the hidden differential genes under the overlap
  if (mode == "shared" && r$differential > 0)
    g <- g + txt(xo, -1.9*lh, sprintf("(+%d\ndiff.)", r$differential), DIFc,
                 2.7 * opt$font_scale, "italic", lineheight=0.85)
  g
}

# --------------------------------------------------------------- run -----------
tab <- bind_rows(lapply(TR, counts_for))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
write_tsv(tab, file.path(opt$outdir, "venn_counts.tsv"))
print(as.data.frame(tab))

for (mode in c("shared", "intersection")) {
  md <- file.path(opt$outdir, mode); dir.create(md, showWarnings=FALSE, recursive=TRUE)
  plots <- list()
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]; p <- draw_venn(r, mode); plots[[i]] <- p
    base <- file.path(md, sprintf("venn_%s_%s", r$id,
                                  gsub("[^A-Za-z0-9]+","_", tolower(r$title))))
    ggsave(paste0(base,".pdf"), p, width=opt$width, height=opt$height, device=cairo_pdf)
    ggsave(paste0(base,".png"), p, width=opt$width, height=opt$height, dpi=200)
  }
  # combined 2x2 montage, tagged A-D
  comb <- wrap_plots(plots, ncol=2) + plot_annotation(tag_levels="A")
  ggsave(file.path(md, "venn_combined.pdf"), comb, width=opt$width*2, height=opt$height*2, device=cairo_pdf)
  ggsave(file.path(md, "venn_combined.png"), comb, width=opt$width*2, height=opt$height*2, dpi=200)
  message("wrote 4 ", mode, " figures + combined to ", md)
}
