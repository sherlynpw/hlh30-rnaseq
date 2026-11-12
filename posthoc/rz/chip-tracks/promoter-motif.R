#!/usr/bin/env Rscript
# Zoomed promoter view: HLH-30 ChIP signal + E-box positions.
#
# Plotted in TSS-relative coordinates (upstream negative) so genes on opposite
# strands line up directly. Two motif classes are shown separately:
#   CACGTG - the strict palindromic HLH-30/TFEB E-box (Lin et al. 2018)
#   CANNTG - the generic bHLH E-box; degenerate, ~1 per 256 bp by chance,
#            so its presence alone carries little information.

suppressPackageStartupMessages({
  library(optparse); library(rtracklayer); library(GenomicRanges); library(Biostrings)
  library(BSgenome.Celegans.UCSC.ce11); library(dplyr); library(ggplot2)
  library(patchwork); library(RSQLite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cache",    type = "character", default = "cache"),
  make_option("--gtf",      type = "character", default = "cache/gtf.gtf"),
  make_option("--gene-db",  type = "character", dest = "gene_db", default = "../../../gene_db.sqlite"),
  make_option("--genes",    type = "character", default = "ztf-28,tts-1"),
  make_option("--upstream", type = "integer",   default = 2000),
  make_option("--downstream", type = "integer", default = 500),
  make_option("--accent",   type = "character", default = "#7B3294"),
  make_option("--outdir",   type = "character", default = "out")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
ACC  <- opt$accent
CHIP <- c("SRX331310", "SRX331312", "SRX331314", "SRX331316")
GEN  <- BSgenome.Celegans.UCSC.ce11

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
  GRanges(d[[1]], IRanges(d[[2]] + 1, d[[3]]))
}

gtf   <- import(opt$gtf)
genes <- gtf[gtf$type == "gene"]
con <- dbConnect(SQLite(), opt$gene_db); on.exit(dbDisconnect(con))
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name FROM annotations")
n2w <- setNames(ann$wormbase_gene, ann$gene_name)
PK  <- lapply(CHIP, rdbed)

NBP  <- 600
want <- trimws(strsplit(opt$genes, ",")[[1]])
XLO  <- -opt$upstream; XHI <- opt$downstream

# genomic position -> TSS-relative (upstream negative), strand aware
rel_of <- function(pos, tss, minus) if (minus) tss - pos else pos - tss

collect <- function(nm) {
  g <- genes[genes$gene_id == n2w[[nm]]]; if (!length(g)) return(NULL); g <- g[1]
  minus <- as.character(strand(g)) == "-"
  tss <- if (minus) end(g) else start(g)
  reg <- if (minus) GRanges(to_ucsc(seqnames(g)), IRanges(tss - XHI, tss - XLO))
         else       GRanges(to_ucsc(seqnames(g)), IRanges(tss + XLO, tss + XHI))

  y  <- rowMeans(sapply(CHIP, function(s) { f <- bwf(s); binned(f, reg, NBP) / genome_mean(f) }))
  gp <- seq(start(reg), end(reg), length.out = NBP)
  sig <- tibble(gene = nm, x = rel_of(gp, tss, minus), y = y) %>% arrange(x)

  sq <- getSeq(GEN, reg)[[1]]
  hit <- function(pat) {
    m <- matchPattern(pat, sq, fixed = (pat == "CACGTG"))
    if (!length(m)) return(numeric(0))
    rel_of(start(reg) + start(m) - 1 + 2.5, tss, minus)   # motif centre
  }
  # Strict palindromic HLH-30/TFEB E-box only (Lin 2018; Alliance of Genome
  # Resources). CANNTG is the generic bHLH box and is uninformative here.
  strict <- hit("CACGTG")
  mot <- if (length(strict)) tibble(gene = nm, x = strict) else tibble(gene = character(), x = numeric())

  pk <- do.call(rbind, lapply(PK, function(p) {
    if (!length(p)) return(NULL)
    ov <- subsetByOverlaps(p, reg); if (!length(ov)) return(NULL)
    data.frame(start = start(ov), end = end(ov))
  }))
  if (!is.null(pk) && nrow(pk)) {
    m <- reduce(GRanges("x", IRanges(pk$start, pk$end)))
    a <- rel_of(start(m), tss, minus); b <- rel_of(end(m), tss, minus)
    pk <- data.frame(start = pmin(a, b), end = pmax(a, b))
  }
  list(gene = nm, sig = sig, mot = mot, pk = pk, gene_end = rel_of(if (minus) start(g) else end(g), tss, minus))
}

dat  <- Filter(Negate(is.null), lapply(want, collect))
YMAX <- ceiling(max(vapply(dat, function(d) max(d$sig$y), numeric(1))) / 5) * 5


panel <- function(d, show_axis, show_leg) {
  n_strict <- nrow(d$mot)
  bar_lo <- -YMAX * 0.16; bar_hi <- -YMAX * 0.04

  p <- ggplot(d$sig, aes(x, y)) +
    geom_area(fill = "grey12", colour = NA) +
    annotate("text", x = XLO + 0.03 * (XHI - XLO), y = YMAX * 0.92, label = "HLH-30",
             hjust = 0, size = 3, fontface = "bold", colour = "grey15") +
    coord_cartesian(xlim = c(XLO, XHI), ylim = c(bar_lo, YMAX * 1.05), expand = FALSE) +
    scale_y_continuous(breaks = c(0, YMAX)) +
    labs(x = NULL, y = if (show_axis) "Normalized read counts" else NULL,
         subtitle = sprintf("%s  -  %d x CACGTG in promoter", d$gene, n_strict)) +
    theme_classic(base_size = 10) +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.text.y  = if (show_axis) element_text(size = 7) else element_blank(),
          axis.ticks.y = if (show_axis) element_line(linewidth = .3) else element_blank(),
          axis.line.y  = if (show_axis) element_line(linewidth = .3) else element_blank(),
          axis.title.y = element_text(size = 8, colour = "grey20"),
          plot.subtitle = element_text(size = 9, face = "italic", colour = "grey15"),
          plot.margin = margin(2, 6, 0, 2))
  if (!is.null(d$pk) && nrow(d$pk))
    p <- p + geom_rect(data = d$pk, inherit.aes = FALSE,
                       aes(xmin = start, xmax = end, ymin = bar_lo, ymax = bar_hi), fill = "black")

  # promoter line + motif lollipops + gene start
  g_end <- min(d$gene_end, XHI)
  m <- ggplot() +
    annotate("segment", x = XLO, xend = 0, y = 0, yend = 0, colour = "grey35", linewidth = 0.7) +
    annotate("text", x = XLO / 2, y = 0.42, label = "promoter", size = 2.6, colour = "grey35") +
    annotate("rect", xmin = 0, xmax = g_end, ymin = -0.22, ymax = 0.22, fill = ACC) +
    annotate("segment", x = 0, xend = 0, y = 0.22, yend = 0.78, colour = ACC, linewidth = 0.5) +
    annotate("segment", x = 0, xend = 0.06 * (XHI - XLO), y = 0.78, yend = 0.78, colour = ACC,
             linewidth = 0.5, arrow = arrow(length = unit(3.5, "pt"), type = "open")) +
    theme_classic(base_size = 9)
  if (nrow(d$mot)) {
    m <- m +
      geom_segment(data = d$mot, aes(x = x, xend = x, y = 0, yend = -0.78),
                   colour = ACC, linewidth = 0.6) +
      geom_point(data = d$mot, aes(x = x, y = -0.86), colour = ACC, size = 2.1) +
      annotate("text", x = XHI, y = -0.86, hjust = 1.05, label = "CACGTG",
               colour = ACC, size = 2.4, fontface = "bold")
  } else {
    m <- m + annotate("text", x = XLO / 2, y = -0.7, label = "no CACGTG",
                      colour = "grey45", size = 2.7, fontface = "italic")
  }
  m <- m +
    coord_cartesian(xlim = c(XLO, XHI), ylim = c(-1.15, 0.95), expand = FALSE) +
    scale_x_continuous(breaks = seq(XLO, XHI, by = 1000)) +
    labs(x = "Position relative to TSS (bp)", y = NULL) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 7), axis.title.x = element_text(size = 8),
          legend.position = if (show_leg) "bottom" else "none",
          legend.text = element_text(size = 7.5), legend.key.height = unit(8, "pt"),
          plot.margin = margin(0, 6, 2, 2))
  p / m + plot_layout(heights = c(5, 2.4))
}

plots <- lapply(seq_along(dat), function(i)
  panel(dat[[i]], show_axis = (i == 1), show_leg = (i == length(dat))))
pF <- wrap_plots(plots, nrow = 1)

ggsave(file.path(opt$outdir, "promoter_motif.pdf"), pF, width = 3.5 * length(plots) + 0.5, height = 3.8)
ggsave(file.path(opt$outdir, "promoter_motif.png"), pF, width = 3.5 * length(plots) + 0.5, height = 3.8, dpi = 220)
message("wrote promoter_motif.{pdf,png}")
