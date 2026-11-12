#!/usr/bin/env Rscript
# Count HLH-30 / DAF-16 / PQM-1 sequence elements in a set of promoters and
# render the counts as an annotated heatmap.
#
# Motif provenance (Lin et al. 2018 Nat Commun 9:4400; Tepper et al. 2013 Cell):
#   E-box   CACGTG / CANNTG  bound by HLH-30 (bHLH). Lin's text says HLH-30
#                            sites are enriched for "diverse E-boxes (CANNTG)";
#                            CACGTG is the canonical palindromic core.
#   DBE     TRTTTAC          DAF-16-Bound Element, bound by DAF-16/FOXO.
#   DAE     TGTAAG/TGATAAG   DAF-16-Associated Element. NOT bound by DAF-16 -
#                            bound by PQM-1. Lin found DAEs enriched even at
#                            HLH-30-only sites, implying PQM-1 co-occupancy at
#                            HLH-30 targets.
#
# `expected` is the count expected by chance for that pattern in a window of
# this length, so degenerate motifs can be read against their own background.
#
# Consensus strings are a blunt instrument: CACGTG misses real variant sites
# and CANNTG admits ~1 hit per 256 bp of background. The JASPAR PWM
# (MA1449.1, hlh-30, C. elegans) weights each position instead, so CACATG and
# CATGTG score as genuine sites while random CANNTG matches do not.

suppressPackageStartupMessages({
  library(optparse); library(rtracklayer); library(GenomicRanges); library(Biostrings)
  library(BSgenome.Celegans.UCSC.ce11); library(ggplot2); library(dplyr); library(RSQLite)
  library(TFBSTools); library(motifmatchr); library(jsonlite)
})

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--gtf",        type = "character", default = "cache/gtf.gtf"),
  make_option("--gene-db",    type = "character", dest = "gene_db", default = "../../../gene_db.sqlite"),
  make_option("--genes",      type = "character", default = "tts-1,cpr-4,ztf-28",
              help = "comma-separated; column order follows this order"),
  make_option("--upstream",   type = "integer",   default = 2000),
  make_option("--downstream", type = "integer",   default = 200),
  make_option("--jaspar",     type = "character", default = "MA1449.1",
              help = "JASPAR matrix ID scanned as a PWM (hlh-30)"),
  make_option("--p-cutoff",   type = "double", default = 1e-4, dest = "p_cutoff",
              help = "PWM match p-value threshold, as used by the pipeline scan"),
  make_option("--cache",      type = "character", default = "cache"),
  make_option("--accent",     type = "character", default = "#7B3294"),
  make_option("--denovo",     type = "character", default = "",
              help = "ChIP-derived HLH-30 motif (denovo-motif.R); scanned as a PWM row. '' = skip."),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size; render at the ORIGINAL canvas."),
  make_option("--outdir",     type = "character", default = "out")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
ACC <- opt$accent

# The element name lives in the facet strip with its factor, so the row labels
# carry only the sequence and nothing is repeated down the column.
MOTIFS <- tibble::tribble(
  ~motif,       ~pattern,     ~tf,        ~element,
  "CACGTG",     "CACGTG",     "HLH-30",   "E-box",
  "TRTTTAC",    "TRTTTAC",    "DAF-16",   "DBE",
  "TTGTTTAC",   "TTGTTTAC",   "DAF-16",   "DBE",
  "TGTAAG",     "TGTAAG",     "PQM-1",    "DAE",
  "TGATAAG",    "TGATAAG",    "PQM-1",    "DAE"
)
STRIP <- function(tf, el) paste0(tf, "\n", el)

to_ucsc <- function(x) { x <- as.character(x); ifelse(x == "MtDNA", "chrM", paste0("chr", x)) }
GEN   <- BSgenome.Celegans.UCSC.ce11
gtf   <- import(opt$gtf); genes <- gtf[gtf$type == "gene"]
con <- dbConnect(SQLite(), opt$gene_db); on.exit(dbDisconnect(con))
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name FROM annotations")
n2w <- setNames(ann$wormbase_gene, ann$gene_name)

IUPAC_N <- c(A=1,C=1,G=1,T=1, R=2,Y=2,S=2,W=2,K=2,M=2, B=3,D=3,H=3,V=3, N=4)
is_palindrome <- function(p) p == as.character(reverseComplement(DNAString(p)))
# expected hits by chance: product of per-position probabilities, over all
# offsets, doubled for non-palindromic patterns (both strands searched)
expected_hits <- function(p, len) {
  pr <- prod(IUPAC_N[strsplit(p, "")[[1]]] / 4)
  n_pos <- len - nchar(p) + 1
  pr * n_pos * if (is_palindrome(p)) 1 else 2
}
count_hits <- function(p, sq) {
  n <- length(matchPattern(p, sq, fixed = FALSE))
  if (!is_palindrome(p))
    n <- n + length(matchPattern(reverseComplement(DNAString(p)), sq, fixed = FALSE))
  n
}

# --- JASPAR PWM, cached on first use -----------------------------------------
jf <- file.path(opt$cache, paste0(opt$jaspar, ".json"))
if (!file.exists(jf)) {
  dir.create(opt$cache, showWarnings = FALSE, recursive = TRUE)
  download.file(sprintf("https://jaspar.elixir.no/api/v1/matrix/%s/", opt$jaspar),
                jf, quiet = TRUE)
}
jj  <- fromJSON(jf)
pfm <- do.call(rbind, lapply(c("A", "C", "G", "T"), function(b) as.integer(jj$pfm[[b]])))
rownames(pfm) <- c("A", "C", "G", "T")
PWM <- toPWM(PFMatrix(ID = opt$jaspar, name = jj$name, profileMatrix = pfm))
# Label the PWM row by what it actually accepts, not by its accession: the
# IUPAC consensus of the 6 bp core (positions with any base >=5% collapsed).
IUPAC2 <- c(AC="M", AG="R", AT="W", CG="S", CT="Y", GT="K")
core_cols <- seq(floor((ncol(pfm) - 6) / 2) + 1, length.out = 6)
CORE_CONS <- paste0(vapply(core_cols, function(i) {
  col <- pfm[, i]; keep <- sort(names(col)[col / sum(col) >= 0.05])
  if (length(keep) == 1) keep else {
    k <- IUPAC2[[paste(keep, collapse = "")]]; if (is.null(k)) "N" else k
  }
}, character(1)), collapse = "")
PWM_LAB <- sprintf("~%s", CORE_CONS)
message("core consensus: ", CORE_CONS)
message("JASPAR ", opt$jaspar, " (", jj$name, "), width ", ncol(pfm))

# --- ChIP-derived HLH-30 motif (STREME on the L4 ChIP summits) ---------------
# File: first line ">CONSENSUS", then "A [ counts ]" rows, as denovo-motif.R
# writes it. Labelled by its consensus, marked ** in the figure.
DENOVO <- NULL
if (nzchar(opt$denovo) && file.exists(opt$denovo)) {
  dl <- readLines(opt$denovo)
  dn_pfm <- do.call(rbind, lapply(dl[-1], function(l)
    as.integer(strsplit(trimws(sub(".*\\[(.*)\\].*", "\\1", l)), "\\s+")[[1]])))
  rownames(dn_pfm) <- substr(dl[-1], 1, 1)
  DENOVO <- list(label = sub("^>", "", dl[1]),
                 pwm = toPWM(PFMatrix(ID = "chip_denovo", name = "HLH-30 ChIP de novo",
                                      profileMatrix = dn_pfm[c("A", "C", "G", "T"), ])))
  message("ChIP-derived motif: ", DENOVO$label)
}

want <- trimws(strsplit(opt$genes, ",")[[1]])
rows <- list()
for (nm in want) {
  g <- genes[genes$gene_id == n2w[[nm]]]
  if (!length(g)) { warning("no GTF entry: ", nm); next }
  g <- g[1]
  minus <- as.character(strand(g)) == "-"
  tss <- if (minus) end(g) else start(g)
  reg <- if (minus) GRanges(to_ucsc(seqnames(g)), IRanges(tss - opt$downstream, tss + opt$upstream))
         else       GRanges(to_ucsc(seqnames(g)), IRanges(tss - opt$upstream, tss + opt$downstream))
  sq <- getSeq(GEN, reg)[[1]]
  if (minus) sq <- reverseComplement(sq)
  for (i in seq_len(nrow(MOTIFS))) {
    rows[[length(rows) + 1]] <- tibble(
      gene = nm, motif = MOTIFS$motif[i],
      tf = STRIP(MOTIFS$tf[i], MOTIFS$element[i]),
      pattern = MOTIFS$pattern[i], strand = as.character(strand(g)),
      count = count_hits(MOTIFS$pattern[i], sq),
      expected = round(expected_hits(MOTIFS$pattern[i], length(sq)), 2))
  }
  # PWM scan: expected hits under the p-value threshold are p * positions * 2
  # strands, which is what motifmatchr's cutoff means operationally.
  hits <- matchMotifs(PWM, DNAStringSet(sq), out = "positions",
                      p.cutoff = opt$p_cutoff)[[1]][[1]]
  # matchMotifs reports one hit per strand and the E-box is palindromic, so a
  # single site appears twice; count distinct positions.
  n_pwm <- length(unique(start(hits)))
  rows[[length(rows) + 1]] <- tibble(
    gene = nm, motif = PWM_LAB, tf = STRIP("HLH-30", "E-box"),
    pattern = sprintf("p<%g", opt$p_cutoff), strand = as.character(strand(g)),
    count = n_pwm,
    expected = round(opt$p_cutoff * length(sq), 2))
  if (!is.null(DENOVO)) {
    dh <- matchMotifs(DENOVO$pwm, DNAStringSet(sq), out = "positions",
                      p.cutoff = opt$p_cutoff)[[1]][[1]]
    # overlapping hits on the two strands are one site
    rows[[length(rows) + 1]] <- tibble(
      gene = nm, motif = DENOVO$label, tf = STRIP("HLH-30", "E-box"),
      pattern = sprintf("p<%g", opt$p_cutoff), strand = as.character(strand(g)),
      count = length(reduce(ranges(dh))),
      expected = round(2 * opt$p_cutoff * length(sq), 2))
  }
}
res <- bind_rows(rows) %>% mutate(obs_exp = round(count / pmax(expected, 1e-9), 2))

win_lab <- sprintf("-%d/+%d bp around the TSS", opt$upstream, opt$downstream)
write_tsv <- utils::write.table
write_tsv(res, file.path(opt$outdir, "motif_counts.tsv"), sep = "\t",
          quote = FALSE, row.names = FALSE)
print(as.data.frame(res %>% select(gene, motif, count, expected, obs_exp)))

# Two clean label columns, guaranteed not to overlap:
#   - the group label (TF name + [element]) lives in the LEFT facet strip,
#   - the motif sequence lives in the y-axis text, right against the tiles.
# The strip is its own gtable column, so it cannot collide with the axis text.
# Each group is exactly two rows, so a centred two-line strip drops the TF name
# onto the top row and [element] onto the bottom row. ggtext's element_markdown
# lets the two lines carry different styles (bold accent vs small grey).
library(ggtext)

DN_LAB <- if (is.null(DENOVO)) character(0) else DENOVO$label
ROW_ORDER <- c(MOTIFS$motif[1], PWM_LAB, DN_LAB, MOTIFS$motif[-1])
grp_of <- setNames(c(rep(STRIP("HLH-30", "E-box"), 2 + length(DN_LAB)),
                     STRIP(MOTIFS$tf[-1], MOTIFS$element[-1])), ROW_ORDER)
GRP_LEVELS <- unique(unname(grp_of[ROW_ORDER]))

# markdown strip label: TF name bold+accent on top, [element] small+grey below
strip_md <- setNames(vapply(GRP_LEVELS, function(g) {
  tf <- sub("\n.*", "", g); el <- sub(".*\n", "", g)
  sprintf("<b style='color:%s'>%s</b><br><span style='font-size:%.1fpt;color:#8A8A8F'>[%s]</span>",
          ACC, tf, 8 * opt$font_scale, el)
}, character(1)), GRP_LEVELS)

res <- res %>%
  mutate(gene  = factor(gene, levels = want),
         motif = factor(motif, levels = rev(ROW_ORDER)),
         grp   = factor(unname(grp_of[as.character(motif)]), levels = GRP_LEVELS))

cmax <- max(res$count)
p <- ggplot(res, aes(gene, motif, fill = count)) +
  geom_tile(colour = "white", linewidth = 1.4) +
  geom_text(aes(label = count,
                colour = ifelse(count > cmax * 0.55, "hi", "lo")),
            size = 3.6 * opt$font_scale, fontface = "bold") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(grp = strip_md)) +
  scale_fill_gradientn(colours = PAL$hlh30_seq,
                       limits = c(0, cmax), name = "sites") +
  scale_colour_manual(values = c(hi = "white", lo = "grey25"), guide = "none") +
  # *  = JASPAR MA1449.1 PWM; ** = motif found de novo in the HLH-30 ChIP peaks
  scale_y_discrete(labels = function(v) ifelse(v == PWM_LAB, paste0(v, "*"),
                                        ifelse(v %in% DN_LAB, paste0(v, "**"), v)),
                   expand = expansion(add = 0.15)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11 * opt$font_scale) +
  theme(
    panel.grid         = element_blank(),
    panel.spacing.y    = unit(6, "pt"),
    axis.text.x        = element_text(face = "bold.italic", size = 10 * opt$font_scale, colour = "grey10"),
    axis.text.y        = element_text(family = "mono", size = 8.5 * opt$font_scale, colour = "grey25", hjust = 1,
                                      margin = margin(r = 2)),
    axis.ticks         = element_blank(),
    strip.placement    = "outside",
    strip.text.y.left  = element_markdown(angle = 0, hjust = 1, lineheight = 1.3,
                                          size = 9 * opt$font_scale, halign = 1),
    legend.key.width   = unit(8, "pt"),
    legend.title       = element_text(size = 8 * opt$font_scale, colour = "grey35", hjust = 0.5),   # 10pt, like the keys
    legend.text        = element_text(size = 8 * opt$font_scale, colour = "grey35"),
    plot.margin        = margin(6, 8, 6, 6)
  )

# one row taller when the ChIP-derived motif is added
mc_h <- 3.0 + 0.4 * length(DN_LAB)
ggsave(file.path(opt$outdir, "motif_counts.pdf"), p, width = 4.1, height = mc_h)
ggsave(file.path(opt$outdir, "motif_counts.png"), p, width = 4.1, height = mc_h, dpi = 220)
message("wrote motif_counts.{tsv,pdf,png}")
