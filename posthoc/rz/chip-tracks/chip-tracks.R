#!/usr/bin/env Rscript
# HLH-30 ChIP-seq figures from ChIP-Atlas ce11 (Oki et al. 2018).
#   Panel A : metagene of HLH-30 signal around TSSs, as fold over matched input
#   Panel E : IGV-style locus tracks - gene model + scale bar on top, one
#             coloured coverage lane per experiment below (cf. Oki Fig 3c)
#
# ChIP-Atlas has exactly four HLH-30 experiments for ce11 (2x L4, 2x late
# embryo), each with its own input control. There is no adult HLH-30 ChIP.

suppressPackageStartupMessages({
  library(optparse); library(rtracklayer); library(GenomicRanges)
  library(dplyr); library(ggplot2); library(patchwork); library(RSQLite)
})

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cache",   type = "character", default = "cache"),
  make_option("--gtf",     type = "character", default = "cache/gtf.gtf"),
  make_option("--gene-db", type = "character", dest = "gene_db", default = "../../../gene_db.sqlite"),
  make_option("--genes",   type = "character", default = "ztf-28,tts-1,cpr-4,lgg-2"),
  make_option("--flank",   type = "integer", default = 3000),
  make_option("--nbin",    type = "integer", default = 60),
  make_option("--max-tss", type = "integer", default = 20000, dest = "max_tss"),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size in panel A; render at the ORIGINAL canvas."),
  make_option("--font-scale-tracks", type = "double", default = NA, dest = "font_scale_tracks",
              help = "Same for panel E, which is placed at a different width. Default: --font-scale."),
  make_option("--outdir",  type = "character", default = "out")
)))
if (is.na(opt$font_scale_tracks)) opt$font_scale_tracks <- opt$font_scale
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ChIP experiment -> its matched input control
EXPS <- tibble::tribble(
  ~lane,                 ~chip,       ~input,      ~stage,
  "HLH-30  L4 rep1",     "SRX331310", "SRX331309", "L4",
  "HLH-30  L4 rep2",     "SRX331312", "SRX331311", "L4",
  "HLH-30  emb rep1",    "SRX331314", "SRX331313", "Late embryo",
  "HLH-30  emb rep2",    "SRX331316", "SRX331315", "Late embryo"
)
LANE_COL <- PAL$lane

to_ucsc <- function(x) { x <- as.character(x); ifelse(x == "MtDNA", "chrM", paste0("chr", x)) }
bwf <- function(s) BigWigFile(file.path(opt$cache, paste0(s, ".bw")))

genome_mean <- function(f) {
  si <- seqinfo(f); gr <- GRanges(seqnames(si), IRanges(1, seqlengths(si)))
  s <- summary(f, which = gr, size = 1, type = "mean", defaultValue = 0)
  v <- vapply(s, function(x) as.numeric(score(x))[1], numeric(1))
  sum(v * seqlengths(si), na.rm = TRUE) / sum(as.numeric(seqlengths(si)))
}
binned <- function(f, gr, n) {
  s <- summary(f, which = gr, size = n, type = "mean", defaultValue = 0)
  matrix(unlist(lapply(s, function(x) as.numeric(score(x)))), nrow = length(gr), byrow = TRUE)
}

message("Reading GTF ...")
gtf   <- import(opt$gtf)
genes <- gtf[gtf$type == "gene"]
pc    <- genes[genes$gene_biotype %in% "protein_coding"]
tss_of <- function(gr) {
  s <- ifelse(as.character(strand(gr)) == "-", end(gr), start(gr))
  GRanges(to_ucsc(seqnames(gr)), IRanges(s, width = 1), strand = strand(gr))
}

# ------------------------------------------------------------------ Panel A
tss <- tss_of(pc)
if (length(tss) > opt$max_tss) tss <- tss[sort(sample(length(tss), opt$max_tss))]
win <- suppressWarnings(trim(resize(tss, width = 2 * opt$flank + 1, fix = "center")))
win <- win[width(win) == 2 * opt$flank + 1]
neg <- as.character(strand(win)) == "-"
message("Metagene over ", length(win), " TSSs (fold over matched input)")

prof <- bind_rows(lapply(seq_len(nrow(EXPS)), function(i) {
  fc <- bwf(EXPS$chip[i]); fi <- bwf(EXPS$input[i])
  mc <- binned(fc, win, opt$nbin) / genome_mean(fc)
  mi <- binned(fi, win, opt$nbin) / genome_mean(fi)
  mc[neg, ] <- mc[neg, ncol(mc):1]; mi[neg, ] <- mi[neg, ncol(mi):1]
  fold <- colMeans(mc, na.rm = TRUE) / pmax(colMeans(mi, na.rm = TRUE), 1e-9)
  # Scale so the flanks (|pos| > 2.5 kb) sit at 1: the genome-mean estimator is
  # not comparable between ChIP and input, so only the shape is trustworthy.
  ppos <- seq(-opt$flank, opt$flank, length.out = opt$nbin)
  fold <- fold / mean(fold[abs(ppos) > 2500], na.rm = TRUE)
  message("  ", EXPS$lane[i], "  peak fold = ", signif(max(fold), 3))
  tibble(lane = EXPS$lane[i], stage = EXPS$stage[i],
         pos = seq(-opt$flank, opt$flank, length.out = opt$nbin), fold = fold)
}))
write.csv(prof, file.path(opt$outdir, "metagene_profile.csv"), row.names = FALSE)

pA <- ggplot(prof, aes(pos, fold, colour = lane)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_line(linewidth = 0.85) +
  # the y label already says HLH-30, so the keys drop it and the legend
  # stays narrow enough to clear the peak
  scale_colour_manual(values = LANE_COL, name = NULL,
                      labels = function(x) sub("^HLH-30\\s+", "", x)) +
  scale_x_continuous(breaks = seq(-3000, 3000, 1500)) +
  scale_y_continuous(breaks = scales::breaks_width(0.5)) +
  labs(x = "Position relative to TSS (bp)",
       y = "HLH-30 fold enrichment over input",
       subtitle = sprintf("ChIP-Atlas ce11, n = %s TSSs", format(length(win), big.mark = ","))) +
  theme_bw(base_size = 11 * opt$font_scale) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", .85), colour = NA),
        legend.key.height = unit(10, "pt"), legend.text = element_text(size = 8.8 * opt$font_scale),
        plot.title = element_text(face = "bold", size = 14 * opt$font_scale),
        plot.subtitle = element_text(size = 9 * opt$font_scale, colour = "grey35"))

# ------------------------------------------------------------------ Panel E
con <- dbConnect(SQLite(), opt$gene_db); on.exit(dbDisconnect(con))
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name FROM annotations")
n2w <- setNames(ann$wormbase_gene, ann$gene_name)

rdbed <- function(s) {
  f <- file.path(opt$cache, paste0(s, ".bed"))
  if (!file.exists(f)) return(GRanges())
  d <- read.table(f, sep = "\t", comment.char = "#", quote = "", fill = TRUE, header = FALSE)
  GRanges(d[[1]], IRanges(d[[2]] + 1, d[[3]]), score = suppressWarnings(as.numeric(d[[5]])))
}
PEAKS <- lapply(EXPS$chip, rdbed); names(PEAKS) <- EXPS$lane

NBP <- 800
locus <- function(gname, show_lanes = TRUE) {
  wb <- n2w[[gname]]; g <- genes[genes$gene_id == wb]
  if (!length(g)) { message("skip (no GTF entry): ", gname); return(NULL) }
  g <- g[1]
  tssg <- tss_of(g)
  reg  <- suppressWarnings(trim(resize(tssg, width = 2 * opt$flank + 1, fix = "center")))
  strand(reg) <- "*"
  xs <- seq(start(reg), end(reg), length.out = NBP)

  # Normalised read counts, as in Lin Fig 5e (input is used only in panel a).
  sig <- bind_rows(lapply(seq_len(nrow(EXPS)), function(i) {
    fc <- bwf(EXPS$chip[i])
    tibble(lane = EXPS$lane[i], x = xs,
           y = as.numeric(binned(fc, reg, NBP)) / genome_mean(fc))
  }))
  sig$lane <- factor(sig$lane, levels = EXPS$lane)

  pk <- bind_rows(lapply(EXPS$lane, function(ln) {
    p <- PEAKS[[ln]]; if (!length(p)) return(NULL)
    ov <- subsetByOverlaps(p, reg); if (!length(ov)) return(NULL)
    tibble(lane = ln, start = start(ov), end = end(ov))
  }))

  # gene model + scale bar on top
  ex <- gtf[gtf$type == "exon" & gtf$gene_id == wb]
  sb_len <- 1000
  sb_x0  <- start(reg) + 0.06 * width(reg)
  arrow_dir <- if (as.character(strand(g)) == "-") -1 else 1
  gm <- ggplot() +
    geom_segment(aes(x = sb_x0, xend = sb_x0 + sb_len, y = 2.15, yend = 2.15),
                 linewidth = 0.5, colour = "grey15") +
    annotate("text", x = sb_x0 + sb_len / 2, y = 2.5, label = "1 kb", size = 2.0 * opt$font_scale_tracks, colour = "grey15") +
    geom_segment(aes(x = start(g), xend = end(g), y = 1, yend = 1),
                 colour = "grey35", linewidth = 0.4) +
    geom_rect(data = as.data.frame(ranges(ex)),
              aes(xmin = start, xmax = end, ymin = 0.55, ymax = 1.45), fill = "grey25") +
    geom_segment(aes(x = start(tssg), xend = start(tssg), y = 1.45, yend = 2.0),
                 colour = "#B2182B", linewidth = 0.5) +
    geom_segment(aes(x = start(tssg), xend = start(tssg) + arrow_dir * 0.05 * width(reg),
                     y = 2.0, yend = 2.0), colour = "#B2182B", linewidth = 0.5,
                 arrow = arrow(length = unit(4, "pt"), type = "open")) +
    coord_cartesian(xlim = c(start(reg), end(reg)), ylim = c(0.3, 2.8), expand = FALSE) +
    labs(subtitle = gname) +
    theme_void(base_size = 9 * opt$font_scale_tracks) +
    theme(plot.subtitle = element_text(face = "bold.italic", hjust = 0.5, size = 6 * opt$font_scale_tracks,
                                       margin = margin(b = 2)))

  trk <- ggplot(sig, aes(x, y, fill = lane)) +
    geom_area(colour = NA) +
    # lane names once, on the first locus, without the repeated TF name
    facet_wrap(~ lane, ncol = 1, strip.position = "left",
               labeller = as_labeller(function(x) sub("^HLH-30\\s+", "", x))) +
    scale_fill_manual(values = LANE_COL, guide = "none") +
    scale_y_continuous(breaks = c(0, floor(YMAX / 5) * 5)) +
    coord_cartesian(xlim = c(start(reg), end(reg)),
                    ylim = c(-YMAX * 0.16, YMAX), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 9 * opt$font_scale_tracks) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(), strip.placement = "outside",
          strip.text.y.left = if (show_lanes)
            element_text(angle = 0, size = 5.5 * opt$font_scale_tracks, hjust = 1, colour = "grey25")
            else element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 5.5 * opt$font_scale_tracks),
          panel.spacing.y = unit(1.5, "pt"))
  if (!is.null(pk) && nrow(pk)) {
    pk$lane <- factor(pk$lane, levels = EXPS$lane)
    trk <- trk + geom_rect(data = pk, inherit.aes = FALSE,
                           aes(xmin = start, xmax = end,
                               ymin = -YMAX * 0.13, ymax = -YMAX * 0.03), fill = "black")
  }
  gm / trk + plot_layout(heights = c(1.5, 6))
}

want <- trimws(strsplit(opt$genes, ",")[[1]])
# Pass 1: collect signal from every locus so all panels share one y-scale -
# with free scales the ztf-28/tts-1 contrast disappears.
YMAX <<- 1
pre <- unlist(lapply(want, function(gn) {
  wb <- n2w[[gn]]; g <- genes[genes$gene_id == wb]
  if (!length(g)) return(NULL)
  tssg <- tss_of(g[1])
  reg <- suppressWarnings(trim(resize(tssg, width = 2 * opt$flank + 1, fix = "center")))
  strand(reg) <- "*"
  unlist(lapply(EXPS$chip, function(s) {
    f <- bwf(s); as.numeric(binned(f, reg, NBP)) / genome_mean(f) }))
}))
YMAX <<- as.numeric(quantile(pre, 0.9995, na.rm = TRUE))
message("shared track y-max = ", signif(YMAX, 3))
panels <- Filter(Negate(is.null), lapply(seq_along(want), function(i) locus(want[i], show_lanes = i == 1)))
pE <- wrap_plots(panels, nrow = 1) +
  plot_annotation(caption = paste("HLH-30 ChIP, normalised read counts, shared y-scale across loci.",
                                  "Black bars = MACS2 peak calls (ChIP-Atlas, q<1e-5).", sep = "\n"),
                  theme = theme(plot.caption = element_text(size = 5.5 * opt$font_scale_tracks, colour = "grey40", hjust = 0)))

ggsave(file.path(opt$outdir, "panelA_metagene.pdf"), pA, width = 5.4, height = 4.0)
ggsave(file.path(opt$outdir, "panelA_metagene.png"), pA, width = 5.4, height = 4.0, dpi = 200)
ggsave(file.path(opt$outdir, "panelE_tracks.pdf"), pE, width = 3.1 * length(panels), height = 4.6)
ggsave(file.path(opt$outdir, "panelE_tracks.png"), pE, width = 3.1 * length(panels), height = 4.6, dpi = 200)
message("wrote panels to ", opt$outdir)
