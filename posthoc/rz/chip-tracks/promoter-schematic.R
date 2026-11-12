#!/usr/bin/env Rscript
# Paper-style promoter schematic (no ChIP data):
#   promoter line + ruler, gene block arrow at the TSS, a rounded box for each
#   predicted HLH-30 motif with the TF name above and coordinate below, and a
#   connector to the magnified sequence with the motif bases highlighted.
#
# Wording is deliberately "predicted HLH-30 motif" - these are sequence-scan
# hits (CACGTG), not sites demonstrated by a binding assay.

suppressPackageStartupMessages({
  library(optparse); library(rtracklayer); library(GenomicRanges); library(Biostrings)
  library(BSgenome.Celegans.UCSC.ce11); library(ggplot2); library(ggforce)
  library(dplyr); library(patchwork); library(RSQLite)
  library(TFBSTools); library(motifmatchr); library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--gtf",        type = "character", default = "cache/gtf.gtf"),
  make_option("--gene-db",    type = "character", dest = "gene_db", default = "../../../gene_db.sqlite"),
  make_option("--genes",      type = "character", default = "tts-1,ztf-28"),
  make_option("--upstream",   type = "integer",   default = 2000),
  make_option("--downstream", type = "integer",   default = 200),
  make_option("--motif-source", type = "character", default = "consensus", dest = "motif_source",
              help = "pwm (JASPAR PFM, default) or consensus (exact string match)"),
  make_option("--motif",      type = "character", default = "CACGTG",
              help = "consensus string, used when --motif-source consensus"),
  make_option("--jaspar",     type = "character", default = "MA1449.1"),
  make_option("--p-cutoff",   type = "double", default = 1e-4, dest = "p_cutoff"),
  make_option("--cache",      type = "character", default = "cache"),
  make_option("--flankbp",    type = "integer",   default = 9),
  make_option("--seqsize",    type = "double",    default = 6.0,
              help = "font size of the magnified sequence"),
  make_option("--seqwidth",   type = "double",    default = 0.60,
              help = "fraction of the panel width the sequence occupies"),
  make_option("--accent",     type = "character", default = "#7B3294"),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size; render at the ORIGINAL canvas."),
  make_option("--outdir",     type = "character", default = "out")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
ACC <- opt$accent; MOT <- opt$motif
ANN <- 3.1 * opt$font_scale    # ruler / caption / coordinate text size
GEN <- BSgenome.Celegans.UCSC.ce11
USE_PWM <- identical(opt$motif_source, "pwm")

# The PWM is 12 bp; only its central 6 bp E-box core is highlighted, so the
# figure reads the same whichever source is used.
CORE <- 6L
if (USE_PWM) {
  jf <- file.path(opt$cache, paste0(opt$jaspar, ".json"))
  if (!file.exists(jf)) {
    dir.create(opt$cache, showWarnings = FALSE, recursive = TRUE)
    download.file(sprintf("https://jaspar.elixir.no/api/v1/matrix/%s/", opt$jaspar), jf, quiet = TRUE)
  }
  jj  <- fromJSON(jf)
  pfm <- do.call(rbind, lapply(c("A","C","G","T"), function(b) as.integer(jj$pfm[[b]])))
  rownames(pfm) <- c("A","C","G","T")
  PWM <- toPWM(PFMatrix(ID = opt$jaspar, name = jj$name, profileMatrix = pfm))
  CORE_OFF <- floor((ncol(pfm) - CORE) / 2)
  CAPTION  <- sprintf("predicted HLH-30 motif (JASPAR %s, p<%g)", opt$jaspar, opt$p_cutoff)
  NONE_LAB <- sprintf("no predicted HLH-30 motif (JASPAR %s, p<%g)", opt$jaspar, opt$p_cutoff)
} else {
  CAPTION  <- sprintf("predicted HLH-30 motif (%s)", MOT)
  NONE_LAB <- sprintf("no predicted HLH-30 motif (%s)", MOT)
}

to_ucsc <- function(x) { x <- as.character(x); ifelse(x == "MtDNA", "chrM", paste0("chr", x)) }
gtf   <- import(opt$gtf); genes <- gtf[gtf$type == "gene"]
con <- dbConnect(SQLite(), opt$gene_db); on.exit(dbDisconnect(con))
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name FROM annotations")
n2w <- setNames(ann$wormbase_gene, ann$gene_name)

XLO <- -opt$upstream; XHI <- opt$downstream
SPAN <- XHI - XLO

collect <- function(nm) {
  g <- genes[genes$gene_id == n2w[[nm]]]; if (!length(g)) return(NULL); g <- g[1]
  minus <- as.character(strand(g)) == "-"
  tss <- if (minus) end(g) else start(g)
  reg <- if (minus) GRanges(to_ucsc(seqnames(g)), IRanges(tss - XHI, tss - XLO))
         else       GRanges(to_ucsc(seqnames(g)), IRanges(tss + XLO, tss + XHI))
  sq  <- getSeq(GEN, reg)[[1]]
  if (minus) sq <- reverseComplement(sq)     # orient 5'->3' along transcription
  # start position of the 6 bp core, whichever source is used
  cs <- if (USE_PWM) {
    h <- matchMotifs(PWM, DNAStringSet(sq), out = "positions", p.cutoff = opt$p_cutoff)[[1]][[1]]
    # one hit per strand for a palindromic motif - keep distinct positions only
    if (length(h)) sort(unique(start(h))) + CORE_OFF else integer(0)
  } else {
    m <- matchPattern(MOT, sq); if (length(m)) start(m) else integer(0)
  }
  pos <- if (length(cs)) XLO + cs - 1 else numeric(0)
  ctx <- if (length(cs)) vapply(seq_along(cs), function(i) {
    a <- max(1, cs[i] - opt$flankbp); b <- min(length(sq), cs[i] + CORE - 1 + opt$flankbp)
    as.character(subseq(sq, a, b))
  }, character(1)) else character(0)
  off <- if (length(cs)) cs - pmax(1, cs - opt$flankbp) else integer(0)
  list(gene = nm, pos = pos, ctx = ctx, off = off)
}

dat <- Filter(Negate(is.null), lapply(trimws(strsplit(opt$genes, ",")[[1]]), collect))

panel <- function(d, tss_below = TRUE) {
  n <- length(d$pos)
  seq_y <- if (n) -1.08 - (seq_len(n) - 1) * 0.80 else numeric(0)
  ybot  <- if (n) min(seq_y) - 0.76 else -1.12
  ytop  <- if (tss_below) 0.54 else 1.08

  p <- ggplot() +
    annotate("segment", x = XLO, xend = 0, y = 0, yend = 0,
             colour = "grey30", linewidth = 1.1) +
    annotate("text", x = XLO - 0.015 * SPAN, y = 0, label = "5'", hjust = 1,
             size = 3 * opt$font_scale, colour = "grey30") +
    annotate("segment", x = seq(XLO, 0, by = 500), xend = seq(XLO, 0, by = 500),
             y = -0.12, yend = -0.26, colour = "grey55", linewidth = 0.35) +
    annotate("text", x = seq(XLO, 0, by = 500), y = -0.45,
             label = format(seq(XLO, 0, by = 500)), size = ANN, colour = "grey45") +
    geom_shape(data = data.frame(
      x = c(0, XHI * 0.55, XHI * 0.55, XHI, XHI * 0.55, XHI * 0.55, 0),
      y = c(-0.26, -0.26, -0.44, 0, 0.44, 0.26, 0.26)),
      aes(x, y), fill = ACC, colour = NA, radius = unit(1, "pt")) +
    annotate("text", x = XHI * 1.08, y = 0, label = d$gene, hjust = 0,
             size = 3.4 * opt$font_scale, fontface = "bold.italic", colour = "grey10") +
    # TSS tag: a white-on-accent pill, under the first row's arrow and over
    # the last row's, so both sit between the two lines and cost no height
    annotate("label", x = XHI * 0.95, y = if (tss_below) -0.84 else 0.78, label = "TSS",
             size = 3.0 * opt$font_scale, fontface = "bold", colour = "white", fill = ACC,
             label.size = 0, label.r = unit(8, "pt"), label.padding = unit(4.5, "pt"))

  if (n) {
    box_w <- 0.011 * SPAN
    # drawn first so the sequence letters render over the dashes
    p <- p + geom_segment(data = data.frame(x = d$pos, ye = seq_y + 0.05),
                          aes(x = x, xend = x, y = -0.30, yend = ye),
                          colour = ACC, linewidth = 0.3, linetype = "22")

    for (i in seq_len(n)) {
      x0 <- d$pos[i]

      ch <- strsplit(d$ctx[i], "")[[1]]; L <- length(ch)
      w  <- opt$seqwidth * SPAN
      xs <- seq(x0 - w / 2, x0 + w / 2, length.out = L)
      xs <- xs + max(0, XLO - min(xs)) - max(0, max(xs) - XHI)
      is_mot <- seq_len(L) >= (d$off[i] + 1) & seq_len(L) <= (d$off[i] + CORE)

      p <- p +
        geom_shape(data = data.frame(x = c(x0 - box_w, x0 + box_w, x0 + box_w, x0 - box_w),
                                     y = c(-0.22, -0.22, 0.22, 0.22)),
                   aes(x, y), fill = ACC, colour = NA, radius = unit(2.2, "pt")) +
        annotate("text", x = xs, y = seq_y[i], label = ch, family = "mono",
                 size = opt$seqsize, fontface = ifelse(is_mot, "bold", "plain"),
                 colour = ifelse(is_mot, ACC, "grey45")) +
        annotate("text", x = min(xs) - 0.028 * SPAN, y = seq_y[i], label = "5'",
                 size = ANN, colour = "grey55", hjust = 1) +
        annotate("text", x = min(xs) - 0.075 * SPAN, y = seq_y[i],
                 label = sprintf("%+d", round(x0)), hjust = 1,
                 size = ANN, colour = ACC) +
        annotate("text", x = max(xs) + 0.028 * SPAN, y = seq_y[i], label = "3'",
                 size = ANN, colour = "grey55", hjust = 0) +
        annotate("segment", x = min(xs[is_mot]) - 0.006 * SPAN,
                 xend = max(xs[is_mot]) + 0.006 * SPAN,
                 y = seq_y[i] - 0.26, yend = seq_y[i] - 0.26,
                 colour = ACC, linewidth = 0.4)

      if (i == n)
        p <- p + annotate("text", x = mean(xs[is_mot]), y = seq_y[i] - 0.52,
                          label = CAPTION,
                          size = ANN, colour = ACC)
    }
  } else {
    p <- p + annotate("text", x = XLO / 2, y = -0.90,
                      label = sprintf("%s in %d bp upstream", NONE_LAB, opt$upstream),
                      size = ANN, colour = "grey45", fontface = "italic")
  }

  p + coord_cartesian(xlim = c(XLO - 0.17 * SPAN, XHI + 0.16 * SPAN),
                      ylim = c(ybot, ytop), expand = FALSE) +
    theme_void()
}

plots <- lapply(seq_along(dat), function(i) panel(dat[[i]], tss_below = i == 1))
hs <- vapply(dat, function(d) max(0.90, 0.62 + 0.78 * length(d$pos)), numeric(1))
pS <- wrap_plots(plots, ncol = 1, heights = hs)

ggsave(file.path(opt$outdir, "promoter_schematic.pdf"), pS, width = 7.6, height = 3.25,
       device = cairo_pdf)  # cairo_pdf resolves the "mono" sequence font like the PNG device
ggsave(file.path(opt$outdir, "promoter_schematic.png"), pS, width = 7.6, height = 3.25, dpi = 220)
message("wrote promoter_schematic.{pdf,png}")
