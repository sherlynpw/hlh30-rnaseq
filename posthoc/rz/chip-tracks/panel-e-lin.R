#!/usr/bin/env Rscript
# Panel e in the layout of Lin et al. 2018 Fig 5e:
#   one HLH-30 track per locus (black filled coverage), MACS2 peak bar beneath,
#   gene models at the BOTTOM (all genes in the window, greedily packed into
#   rows), "<gene> promoter" caption underneath, and a 1 kb scale bar.
# Single shared y-axis ("Normalized coverage", coverage over its genome mean) across all loci.

suppressPackageStartupMessages({
  library(optparse); library(rtracklayer); library(GenomicRanges)
  library(dplyr); library(ggplot2); library(patchwork); library(RSQLite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cache",   type = "character", default = "cache"),
  make_option("--gtf",     type = "character", default = "cache/gtf.gtf"),
  make_option("--gene-db", type = "character", dest = "gene_db", default = "../../../gene_db.sqlite"),
  make_option("--genes",   type = "character", default = "ztf-28,tts-1,cpr-4"),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size; render at the ORIGINAL canvas."),
  make_option("--stages",  type = "character", default = "all",
              help = "'all' = mean of 4 experiments; 'L4' = mean of the two L4 reps"),
  make_option("--flank",   type = "integer",   default = 2500),
  make_option("--accent",  type = "character", default = "#7B3294",
              help = "gene-model colour"),
  make_option("--outdir",  type = "character", default = "out")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

ALL  <- c("SRX331310", "SRX331312", "SRX331314", "SRX331316")
CHIP <- if (identical(opt$stages, "L4")) ALL[1:2] else ALL
ACC  <- opt$accent

to_ucsc <- function(x) { x <- as.character(x); ifelse(x == "MtDNA", "chrM", paste0("chr", x)) }
bwf <- function(s) BigWigFile(file.path(opt$cache, paste0(s, ".bw")))
genome_mean <- function(f) {
  si <- seqinfo(f); gr <- GRanges(seqnames(si), IRanges(1, seqlengths(si)))
  s <- summary(f, which = gr, size = 1, type = "mean", defaultValue = 0)
  v <- vapply(s, function(x) as.numeric(score(x))[1], numeric(1))
  sum(v * seqlengths(si), na.rm = TRUE) / sum(as.numeric(seqlengths(si)))
}
binned <- function(f, gr, n) as.numeric(score(summary(f, which = gr, size = n,
                                                      type = "mean", defaultValue = 0)[[1]]))
rdbed <- function(s) {
  f <- file.path(opt$cache, paste0(s, ".bed")); if (!file.exists(f)) return(GRanges())
  d <- read.table(f, sep = "\t", comment.char = "#", quote = "", fill = TRUE, header = FALSE)
  GRanges(d[[1]], IRanges(d[[2]] + 1, d[[3]]), score = suppressWarnings(as.numeric(d[[5]])))
}

gtf   <- import(opt$gtf)
genes <- gtf[gtf$type == "gene"]
con <- dbConnect(SQLite(), opt$gene_db); on.exit(dbDisconnect(con))
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name FROM annotations")
n2w <- setNames(ann$wormbase_gene, ann$gene_name)
w2n <- setNames(ann$gene_name, ann$wormbase_gene)
PK  <- lapply(CHIP, rdbed)

NBP  <- 700
want <- trimws(strsplit(opt$genes, ",")[[1]])

dat <- lapply(want, function(gn) {
  g <- genes[genes$gene_id == n2w[[gn]]]; if (!length(g)) return(NULL); g <- g[1]
  t0  <- ifelse(as.character(strand(g)) == "-", end(g), start(g))
  reg <- GRanges(to_ucsc(seqnames(g)), IRanges(t0 - opt$flank, t0 + opt$flank))
  y <- rowMeans(sapply(CHIP, function(s) { f <- bwf(s); binned(f, reg, NBP) / genome_mean(f) }))
  pk <- do.call(rbind, lapply(PK, function(p) {
    if (!length(p)) return(NULL)
    ov <- subsetByOverlaps(p, reg); if (!length(ov)) return(NULL)
    data.frame(start = start(ov), end = end(ov))
  }))
  if (!is.null(pk) && nrow(pk)) {
    m <- reduce(GRanges("x", IRanges(pk$start, pk$end)))
    pk <- data.frame(start = start(m), end = end(m))
  }
  # every gene in the window, not just the target
  regE <- GRanges(as.character(seqnames(g)), ranges(reg))
  nb <- subsetByOverlaps(genes, regE)
  list(name = gn, g = g, tss = t0, reg = reg, nb = nb, pk = pk,
       x = seq(start(reg), end(reg), length.out = NBP), y = y)
})
dat  <- Filter(Negate(is.null), dat)
YMAX <- ceiling(max(vapply(dat, function(d) max(d$y), numeric(1))) / 5) * 5
message("shared y-max = ", YMAX)

# greedy packing so overlapping genes get their own row
pack_rows <- function(st, en, pad) {
  ord <- order(st); row <- integer(length(st)); ends <- numeric(0)
  for (i in ord) {
    r <- which(ends < st[i] - pad)[1]
    if (is.na(r)) { ends <- c(ends, en[i]); r <- length(ends) } else ends[r] <- en[i]
    row[i] <- r
  }
  row
}

# Every locus reserves the same number of gene-model rows, so the track panels
# are all the same height and their x-axes line up across loci.
NROWMAX <- max(vapply(dat, function(d)
  max(pack_rows(start(d$nb), end(d$nb), pad = 0.16 * width(d$reg))), numeric(1)))
message("gene-model rows reserved per panel: ", NROWMAX)

panel <- function(d, show_axis) {
  bar_lo <- -YMAX * 0.17; bar_hi <- -YMAX * 0.06
  sb_len <- 1000
  sb_x1  <- end(d$reg) - 0.04 * width(d$reg); sb_x0 <- sb_x1 - sb_len

  p <- ggplot(data.frame(x = d$x, y = d$y), aes(x, y)) +
    geom_area(fill = "grey12", colour = NA) +
    annotate("segment", x = sb_x0, xend = sb_x1, y = YMAX * 0.93, yend = YMAX * 0.93,
             linewidth = 0.6, colour = "grey10") +
    annotate("text", x = (sb_x0 + sb_x1) / 2, y = YMAX * 0.99, label = "1 kb",
             size = 2.7 * opt$font_scale, colour = "grey10", vjust = 0.6) +
    annotate("text", x = start(d$reg) + 0.03 * width(d$reg), y = YMAX * 0.93,
             label = "HLH-30", hjust = 0, size = 3 * opt$font_scale, fontface = "bold", colour = "grey15") +
    coord_cartesian(xlim = c(start(d$reg), end(d$reg)),
                    ylim = c(bar_lo, YMAX * 1.04), expand = FALSE) +
    scale_y_continuous(breaks = c(0, YMAX)) +
    labs(x = NULL, y = if (show_axis) "Normalized coverage" else NULL) +
    theme_classic(base_size = 10 * opt$font_scale) +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y  = if (show_axis) element_text(size = 7 * opt$font_scale) else element_blank(),
          axis.ticks.y = if (show_axis) element_line(linewidth = .3) else element_blank(),
          axis.line.y  = if (show_axis) element_line(linewidth = .3) else element_blank(),
          axis.title.y = element_text(size = 8 * opt$font_scale, colour = "grey20"),
          plot.margin = margin(2, 4, 0, 2))
  if (!is.null(d$pk) && nrow(d$pk))
    p <- p + geom_rect(data = d$pk, inherit.aes = FALSE,
                       aes(xmin = start, xmax = end, ymin = bar_lo, ymax = bar_hi),
                       fill = "black")

  nb  <- d$nb
  lab <- ifelse(is.na(w2n[nb$gene_id]) | w2n[nb$gene_id] == "", nb$gene_id, w2n[nb$gene_id])
  gstr <- as.character(strand(nb))
  gdf <- data.frame(start = start(nb), end = end(nb), lab = lab,
                    strand = gstr,
                    tss_x = ifelse(gstr == "-", end(nb), start(nb)),
                    dir   = ifelse(gstr == "-", -1, 1),
                    is_target = nb$gene_id == d$g$gene_id,
                    stringsAsFactors = FALSE)
  gdf$row <- pack_rows(gdf$start, gdf$end, pad = 0.16 * width(d$reg))
  ROW_H <- 1.4                      # row pitch in model units; labels need the air at paper type size
  gdf$yc  <- -gdf$row * ROW_H
  # Genes that extend past the window have their midpoint off-screen; clamp the
  # label so it stays visible.
  gdf$col  <- ifelse(gdf$is_target, ACC, "grey62")
  gdf$lcol <- ifelse(gdf$is_target, "grey10", "grey55")
  gdf$lx <- pmin(pmax((gdf$start + gdf$end) / 2,
                      start(d$reg) + 0.10 * width(d$reg)),
                 end(d$reg) - 0.10 * width(d$reg))
  exdf <- do.call(rbind, lapply(seq_len(nrow(gdf)), function(i) {
    ex <- gtf[gtf$type == "exon" & gtf$gene_id == nb$gene_id[i]]
    if (!length(ex)) return(NULL)
    data.frame(start = start(ex), end = end(ex), yc = gdf$yc[i], col = gdf$col[i])
  }))
  nrow_max <- NROWMAX

  gm <- ggplot() +
    geom_segment(data = gdf, aes(x = start, xend = end, y = yc, yend = yc, colour = col),
                 linewidth = 0.4) +
    geom_rect(data = exdf, aes(xmin = start, xmax = end, ymin = yc - 0.26, ymax = yc + 0.26,
                               fill = col)) +
    geom_text(data = gdf, aes(x = lx, y = yc - 0.62, label = lab,
                              fontface = ifelse(is_target, "bold.italic", "italic"),
                              colour = lcol), size = 2.1 * opt$font_scale) +
    # small arrow at every gene's TSS, pointing in the direction of transcription
    geom_segment(data = gdf, aes(x = tss_x, xend = tss_x, y = yc + 0.26, yend = yc + 0.52,
                                 colour = col), linewidth = 0.32) +
    geom_segment(data = gdf,
                 aes(x = tss_x, xend = tss_x + dir * 0.055 * width(d$reg),
                     y = yc + 0.52, yend = yc + 0.52, colour = col),
                 linewidth = 0.32,
                 arrow = arrow(length = unit(2.6, "pt"), type = "open")) +
    scale_colour_identity() + scale_fill_identity() +
    coord_cartesian(xlim = c(start(d$reg), end(d$reg)),
                    ylim = c(-nrow_max * ROW_H - 0.95, 0.1), expand = FALSE) +
    labs(caption = paste0(d$name, " promoter")) +
    theme_void(base_size = 9 * opt$font_scale) +
    theme(plot.caption = element_text(hjust = 0.5, size = 9.5 * opt$font_scale, face = "italic",
                                      colour = "grey10", margin = margin(t = 3)),
          plot.margin = margin(0, 4, 2, 2))
  p / gm + plot_layout(heights = c(6, 1.1 + 0.7 * nrow_max * ROW_H))
}

plots <- lapply(seq_along(dat), function(i) panel(dat[[i]], show_axis = (i == 1)))
# no panel letter - the position in the final figure is not settled yet
pE <- wrap_plots(plots, nrow = 1)

# taller with the type: the gene-model rows and their labels are in text units
ggsave(file.path(opt$outdir, "panelE_lin_style.pdf"), pE, width = 2.6 * length(plots) + 0.4, height = 3.6 * (0.6 + 0.4 * opt$font_scale))
ggsave(file.path(opt$outdir, "panelE_lin_style.png"), pE, width = 2.6 * length(plots) + 0.4, height = 3.6 * (0.6 + 0.4 * opt$font_scale), dpi = 220)
message("wrote panelE_lin_style.{pdf,png}")
