#!/usr/bin/env Rscript
# What fraction of each gene list carries an HLH-30 motif in the promoter
# -2000/+200 of the TSS, tested against the expressed-gene background with
# Fisher's exact test. Four motif definitions are scanned, plus the union of
# the two HLH-30 motifs:
#   any     - a site of hlh30 or chip
#   hlh30   - canonical HLH-30 E-box, JASPAR MA1449.1 PWM (TCACGTGA)
#   chip    - the motif found de novo in the HLH-30 ChIP peaks (STREME,
#             chip-tracks/denovo-motif.R, CACGTGRT)
#   dbe     - a DAF-16-binding element (DBE, TRTTTAC, either strand). HLH-30
#             can bind DNA as a complex with DAF-16 (Lin et al. 2018), so a
#             partner site alone may hold HLH-30 without an E-box
#   dae     - a DAF-16-associated element (DAE, CTTATCA, either strand), the
#             PQM-1 site (Tepper et al. 2013)
# The scan also keeps cob_dbe / cob_dae: hlh30 sites with that element within
# --cobound-bp (not plotted).
#
# This is the motif-based twin of chip-promoter.R (which uses experimental
# ChIP-Atlas peaks). Motif presence is a genomic property, so the promoter is
# scanned once for the whole expressed universe and cached; every gene list is
# then a lookup + a 2x2 test. Gene lists come from three families:
#   deg      - DESeq2 DEGs per timepoint x {all, up, down}
#   compare  - compare_de category lists (exclusive / shared / differential)
#   go       - genes appearing in significantly enriched GO terms of a DEG list
#
# Emits motif_targets_summary.tsv; render-targets.R draws the figures.

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(tidyr); library(RSQLite)
  library(rtracklayer); library(GenomicRanges); library(Biostrings)
  library(BSgenome.Celegans.UCSC.ce11); library(TFBSTools); library(motifmatchr); library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run",        type = "character", default = "hlh30-all-ruv"),
  make_option("--gtf",        type = "character", default = "cache/gtf.gtf"),
  make_option("--gene-db",    type = "character", dest = "gene_db", default = "../../../gene_db.sqlite"),
  make_option("--de-dir",     type = "character", dest = "de_dir",
              default = "../chip-promoter/cache/hlh30-all-ruv/de",
              help = "DESeq2 <contrast>.csv dir (defines the expressed background)"),
  make_option("--deg-dir",    type = "character", dest = "deg_dir",
              default = "../chip-promoter/cache/hlh30-all-ruv/deg",
              help = "dir of <contrast>_genes_{de,up,down}.txt"),
  make_option("--compare-dir", type = "character", dest = "compare_dir",
              default = "../compare-de/out/hlh30-all-ruv/annotated",
              help = "dir of compare_de annotated *.tsv (gene_id + direction)"),
  make_option("--go-dir",     type = "character", dest = "go_dir", default = "cache/go",
              help = "optional dir of GO CSVs, mirroring <list>/go_<ONT>.csv"),
  make_option("--contrasts",  type = "character",
              default = "tm1978_L1_vs_WT_L1,tm1978_L4_vs_WT_L4,tm1978_Day3_vs_WT_Day3,tm1978_Day6_vs_WT_Day6"),
  make_option("--jaspar",     type = "character", default = "MA1449.1"),
  make_option("--chip-motif", type = "character", dest = "chip_motif",
              default = "../chip-tracks/out/denovo/chip_motif.jaspar"),
  make_option("--dbe",        type = "character", default = "TRTTTAC"),
  make_option("--dae",        type = "character", default = "CTTATCA"),
  make_option("--cobound-bp", type = "integer", default = 500, dest = "cobound_bp"),
  make_option("--sets",       type = "character", default = "any,hlh30,chip",
              help = "motif sets to test. dbe and dae alone sit in 36-52%% of all promoters, so they cannot mark targets"),
  make_option("--p-cutoff",   type = "double", default = 1e-4, dest = "p_cutoff"),
  make_option("--upstream",   type = "integer", default = 2000),
  make_option("--downstream", type = "integer", default = 200),
  make_option("--go-padj",    type = "double", default = 0.05, dest = "go_padj"),
  make_option("--cache",      type = "character", default = "cache"),
  make_option("--outdir",     type = "character", default = "out")
)))
contrasts <- trimws(strsplit(opt$contrasts, ",")[[1]])
outdir <- file.path(opt$outdir, opt$run); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$cache, recursive = TRUE, showWarnings = FALSE)

# --- annotation maps ---------------------------------------------------------
con <- dbConnect(SQLite(), opt$gene_db); on.exit(dbDisconnect(con))
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name, sequence_name FROM annotations")
sym_of <- ifelse(!is.na(ann$gene_name) & nzchar(ann$gene_name), ann$gene_name, ann$sequence_name)
# symbol -> WBGene, using both gene_name and sequence_name so GO symbols resolve
sym2wb <- c(setNames(ann$wormbase_gene, ann$gene_name),
            setNames(ann$wormbase_gene, ann$sequence_name))
sym2wb <- sym2wb[!is.na(names(sym2wb)) & nzchar(names(sym2wb))]

# --- expressed universe = union of tested (padj not NA) genes ----------------
universe <- unique(unlist(lapply(contrasts, function(ct) {
  de <- read_csv(file.path(opt$de_dir, paste0(ct, ".csv")), show_col_types = FALSE)
  de$gene_id[!is.na(de$padj)]
})))
message("expressed universe: ", length(universe), " genes")

# --- promoter scan over the universe (cached) ---------------------------------
to_ucsc <- function(x) { x <- as.character(x); ifelse(x == "MtDNA", "chrM", paste0("chr", x)) }
MOTIFS <- c(any = "Any HLH-30 motif", hlh30 = "HLH-30 E-box (JASPAR)", chip = "ChIP-seq motif",
            dbe = "DAF-16 site (DBE)", dae = "PQM-1 site (DAE)")
MOTIFS <- MOTIFS[trimws(strsplit(opt$sets, ",")[[1]])]
scan_cache <- file.path(opt$cache, sprintf("motif_scan_v2_%s_u%d_d%d_p%g_%s_%s%d.tsv",
                                           opt$jaspar, opt$upstream, opt$downstream,
                                           opt$p_cutoff, opt$dbe, opt$dae, opt$cobound_bp))
if (file.exists(scan_cache)) {
  hm <- read_tsv(scan_cache, show_col_types = FALSE)
  message("loaded cached motif scan: ", scan_cache)
} else {
  jf <- file.path(opt$cache, paste0(opt$jaspar, ".json"))
  if (!file.exists(jf))
    download.file(sprintf("https://jaspar.elixir.no/api/v1/matrix/%s/", opt$jaspar), jf, quiet = TRUE)
  jj  <- fromJSON(jf)
  pfm <- do.call(rbind, lapply(c("A","C","G","T"), function(b) as.integer(jj$pfm[[b]])))
  rownames(pfm) <- c("A","C","G","T")
  PWM <- toPWM(PFMatrix(ID = opt$jaspar, name = jj$name, profileMatrix = pfm))
  # de novo ChIP motif, JASPAR-style count matrix written by denovo-motif.R
  dl <- readLines(opt$chip_motif)
  cm <- do.call(rbind, lapply(dl[-1], function(l)
    as.integer(strsplit(trimws(sub(".*\\[(.*)\\].*", "\\1", l)), "\\s+")[[1]])))
  rownames(cm) <- substr(dl[-1], 1, 1)
  CHIP <- toPWM(PFMatrix(ID = "chip", profileMatrix = cm[c("A","C","G","T"), ]))

  gtf   <- import(opt$gtf); genes <- gtf[gtf$type == "gene"]
  gi    <- genes$gene_id
  keep  <- intersect(universe, gi)
  g     <- genes[match(keep, gi)]
  minus <- as.character(strand(g)) == "-"
  tss   <- ifelse(minus, end(g), start(g))
  reg   <- GRanges(to_ucsc(seqnames(g)),
                   IRanges(ifelse(minus, tss - opt$downstream, tss - opt$upstream),
                           ifelse(minus, tss + opt$upstream,   tss + opt$downstream)))
  # clamp windows that run off a chromosome end before fetching
  gen_si <- seqinfo(BSgenome.Celegans.UCSC.ce11)
  seqlevels(reg) <- seqlevels(gen_si)
  seqinfo(reg) <- gen_si
  start(reg) <- pmax(start(reg), 1L)
  reg   <- trim(reg)
  message("scanning ", length(reg), " promoters ...")
  sq  <- getSeq(BSgenome.Celegans.UCSC.ce11, reg)
  sq[minus] <- reverseComplement(sq[minus])
  hlh  <- matchMotifs(PWM,  sq, out = "positions", p.cutoff = opt$p_cutoff)[[1]]
  chip <- matchMotifs(CHIP, sq, out = "positions", p.cutoff = opt$p_cutoff)[[1]]
  # DBE / DAE starts on both strands; the PWM scans above cover both already
  sites <- function(pat) Map(function(f, r) c(start(f), start(r)),
                             vmatchPattern(pat, sq, fixed = FALSE),
                             vmatchPattern(as.character(reverseComplement(DNAString(pat))), sq, fixed = FALSE))
  # hlh30 sites with an element start within cobound_bp
  near <- function(el) mapply(function(h, d) {
    if (!length(h) || !length(d)) return(0L)
    sum(vapply(start(h), function(x) any(abs(d - x) <= opt$cobound_bp), logical(1)))
  }, as.list(hlh), el)
  dbe <- sites(opt$dbe); dae <- sites(opt$dae)
  hm <- tibble(gene_id = keep, hlh30 = lengths(hlh), chip = lengths(chip),
               dbe = lengths(dbe), cob_dbe = near(dbe),
               dae = lengths(dae), cob_dae = near(dae))
  write_tsv(hm, scan_cache)
}
uni <- intersect(universe, hm$gene_id)
N <- length(uni)
hm$any <- as.integer(hm$hlh30 > 0 | hm$chip > 0)
has <- lapply(names(MOTIFS), function(m) setNames(as.integer(hm[[m]] > 0), hm$gene_id))
names(has) <- names(MOTIFS)
bg <- vapply(has, function(h) 100 * sum(h[uni]) / N, numeric(1))
for (m in names(MOTIFS))
  message(sprintf("background %-8s %d/%d = %.1f%%", m, sum(has[[m]][uni]), N, bg[[m]]))

# --- assemble gene lists -----------------------------------------------------
stage_of <- function(c) c(L1="L1", L4="Day 0", Day3="Day 3", Day6="Day 6")[
  sub("^tm1978_", "", sub("_vs_WT_.*$", "", c))]
lists <- list()   # each: list(id, family, group, direction, wb)
add <- function(id, family, group, direction, wb) {
  wb <- intersect(unique(wb[nzchar(wb)]), uni)
  if (length(wb) >= 3) lists[[length(lists) + 1]] <<-
    list(id = id, family = family, group = group, direction = direction, wb = wb)
}

# deg family
dir_lab <- c(de = "all", up = "up", down = "down")
for (ct in contrasts) for (k in names(dir_lab)) {
  f <- file.path(opt$deg_dir, sprintf("%s_genes_%s.txt", ct, k))
  if (file.exists(f)) add(sprintf("%s | %s", stage_of(ct), dir_lab[k]),
                          "deg", stage_of(ct), dir_lab[k], readLines(f))
}

# compare_de family: split each annotated table by its `direction` column
if (dir.exists(opt$compare_dir)) {
  for (f in list.files(opt$compare_dir, pattern = "\\.tsv$", full.names = TRUE)) {
    base <- sub("\\.tsv$", "", basename(f))
    if (base %in% c("all_categories")) next
    d <- suppressWarnings(read_tsv(f, show_col_types = FALSE))
    if (!all(c("gene_id","direction") %in% names(d))) next
    grp <- sub("_(up|down)$", "", base)
    for (dr in c("up","down")) {
      sub_wb <- d$gene_id[d$direction == dr]
      if (length(sub_wb)) add(sprintf("%s | %s", grp, dr), "compare", grp, dr, sub_wb)
    }
  }
}

# go family: union of genes in enriched terms (p.adjust < go_padj) of a list's GO CSV
if (dir.exists(opt$go_dir)) {
  for (csv in list.files(opt$go_dir, pattern = "go_.*\\.csv$", recursive = TRUE, full.names = TRUE)) {
    go <- suppressWarnings(read_csv(csv, show_col_types = FALSE))
    if (!all(c("geneID","p.adjust") %in% names(go))) next
    sig <- go %>% filter(p.adjust < opt$go_padj)
    if (!nrow(sig)) next
    syms <- unique(unlist(strsplit(paste(sig$geneID, collapse = "/"), "/")))
    wb <- unique(na.omit(sym2wb[syms]))
    # id from the parent list dir, direction inferred from its name
    listname <- basename(dirname(csv))
    ont <- sub("^go_", "", sub("\\.csv$", "", basename(csv)))
    dr <- if (grepl("_up$|_up_|exclusive_up", listname)) "up"
          else if (grepl("_down$|_down_|exclusive_down", listname)) "down" else "all"
    add(sprintf("%s %s | %s", listname, ont, dr), "go", paste(listname, ont), dr, wb)
  }
}
message("assembled ", length(lists), " gene lists")

# --- per-list Fisher enrichment, per motif ---------------------------------
# list vs the rest of the expressed universe. BH runs within each list family
# (the DEG panel's 12 lists x 3 motif sets are one family of 36 tests).
rows <- lapply(lists, function(L) bind_rows(lapply(names(MOTIFS), function(m) {
  h <- has[[m]]; N_hit <- sum(h[uni])
  n <- length(L$wb); a <- sum(h[L$wb])
  b <- n - a; c <- N_hit - a; d <- (N - n) - c
  ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2)),
                 error = function(e) list(p.value = NA_real_, estimate = NA_real_,
                                          conf.int = c(NA_real_, NA_real_)))
  tibble(id = L$id, family = L$family, group = L$group, direction = L$direction,
         motif = m, motif_label = MOTIFS[[m]],
         n = n, n_target = a, pct = 100 * a / n,
         bg_pct = bg[[m]], fold = (a / n) / (N_hit / N),
         odds_ratio = unname(ft$estimate),
         or_lo = ft$conf.int[1], or_hi = ft$conf.int[2], p = ft$p.value)
})))
res <- bind_rows(rows) %>% group_by(family) %>%
  mutate(padj = p.adjust(p, method = "BH")) %>% ungroup() %>%
  mutate(star = ifelse(is.na(padj), "", ifelse(padj < 0.001, "***",
         ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", "")))))

write_tsv(res, file.path(outdir, "motif_targets_summary.tsv"))
print(as.data.frame(res %>% select(id, motif, n, n_target, pct, bg_pct, fold, padj, star)),
      digits = 3, row.names = FALSE)
message("\nwrote ", file.path(outdir, "motif_targets_summary.tsv"))
