#!/usr/bin/env Rscript
# Unbiased motif enrichment in promoters of hlh-30(tm1978) vs WT DEGs, the way
# Chang et al. 2026 (PLoS Biol 24:e3003639) did it: MEME Suite SEA (5.5.7) on
# promoters of significant genes against promoters of non-significant
# (expressed) genes, every motif in JASPAR 2022 nematodes (non-redundant) and
# CIS-BP 2.0 C. elegans, windows 500 and 2,000 bp upstream of the TSS,
# significant at SEA q < 0.05.
#
# One test per stage x direction (up / down in the mutant) x window, plus the
# union over stages. DEGs: padj < 0.05 and |log2FC| >= 1, as in the volcanoes.
#
#   Rscript sea-dge.R              # -> out/sea/{fasta,sea runs,sea_all.tsv}

suppressPackageStartupMessages({
  library(optparse); library(rtracklayer); library(GenomicRanges); library(Biostrings)
  library(BSgenome.Celegans.UCSC.ce11); library(dplyr)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--deg-dir", type = "character", default = "../volcano/cache/hlh30-all-ruv", dest = "deg_dir"),
  make_option("--gtf",     type = "character", default = "../chip-tracks/cache/gtf.gtf"),
  make_option("--padj",    type = "double", default = 0.05),
  make_option("--lfc",     type = "double", default = 1),
  make_option("--windows", type = "character", default = "500,2000"),
  make_option("--image",   type = "character", default = "memesuite/memesuite:5.5.7"),
  make_option("--outdir",  type = "character", default = "out/sea")
)))
dir.create(file.path(opt$outdir, "fasta"), recursive = TRUE, showWarnings = FALSE)
STAGES <- c(L1 = "L1", D0 = "L4", D3 = "Day3", D6 = "Day6")
DB <- "/opt/meme/share/meme-5.5.7/db/motif_databases"
DBS <- c(file.path(DB, "JASPAR", "JASPAR2022_CORE_nematodes_non-redundant_v2.meme"),
         file.path(DB, "CIS-BP_2.00", "Caenorhabditis_elegans.meme"))

gtf <- import(opt$gtf); genes <- gtf[gtf$type == "gene"]
genes <- genes[!duplicated(genes$gene_id)]
names(genes) <- genes$gene_id
to_ucsc <- function(x) { x <- as.character(x); ifelse(x == "MtDNA", "chrM", paste0("chr", x)) }
prom_seq <- function(ids, up) {
  g <- genes[ids]
  minus <- as.character(strand(g)) == "-"; tss <- ifelse(minus, end(g), start(g))
  gr <- GRanges(to_ucsc(seqnames(g)),
                IRanges(ifelse(minus, tss + 1, tss - up), ifelse(minus, tss + up, tss - 1)),
                strand = strand(g))
  names(gr) <- ids
  gr <- gr[as.character(seqnames(gr)) %in% seqnames(BSgenome.Celegans.UCSC.ce11)]
  seqlevels(gr) <- seqlevels(BSgenome.Celegans.UCSC.ce11); seqinfo(gr) <- seqinfo(BSgenome.Celegans.UCSC.ce11)
  gr <- trim(gr)
  s <- getSeq(BSgenome.Celegans.UCSC.ce11, gr)          # strand-aware: 5'->3' of the gene
  names(s) <- names(gr); s
}

sets <- list(); ctrl <- list()
for (st in names(STAGES)) {
  f <- file.path(opt$deg_dir, sprintf("tm1978_%s_vs_WT_%s.csv", STAGES[[st]], STAGES[[st]]))
  d <- read.csv(f, stringsAsFactors = FALSE)
  d <- d[!is.na(d$padj) & d$gene_id %in% names(genes), ]
  deg <- d$padj < opt$padj & abs(d$log2FoldChange) >= opt$lfc
  sets[[paste0(st, "_up")]]   <- d$gene_id[deg & d$log2FoldChange > 0]
  sets[[paste0(st, "_down")]] <- d$gene_id[deg & d$log2FoldChange < 0]
  ctrl[[st]] <- d$gene_id[!is.na(d$padj) & d$padj >= opt$padj]
}
# union over stages; control = expressed at every stage and never a DEG
sets$any_up   <- unique(unlist(sets[grep("_up$", names(sets))]))
sets$any_down <- unique(unlist(sets[grep("_down$", names(sets))]))
ctrl$any <- setdiff(Reduce(intersect, ctrl[names(STAGES)]), c(sets$any_up, sets$any_down))
n <- vapply(sets, length, 1L); print(n)

runs <- list()
for (up in as.integer(strsplit(opt$windows, ",")[[1]])) {
  for (nm in names(sets)) {
    st <- sub("_.*", "", nm)
    if (length(sets[[nm]]) < 10) { message("skip ", nm, " (", length(sets[[nm]]), " genes)"); next }
    tag <- sprintf("%s_u%d", nm, up)
    fp <- file.path(opt$outdir, "fasta", paste0(tag, ".fa"))
    fn <- file.path(opt$outdir, "fasta", sprintf("ctrl_%s_u%d.fa", st, up))
    writeXStringSet(prom_seq(sets[[nm]], up), fp)
    if (!file.exists(fn)) writeXStringSet(prom_seq(ctrl[[st]], up), fn)
    cmd <- sprintf(paste("docker run --rm --user %s -v %s:/work -w /work %s sea --p fasta/%s --n fasta/%s",
                         "%s --oc sea_%s --verbosity 1"),
                   paste(system("id -u", intern = TRUE), system("id -g", intern = TRUE), sep = ":"),
                   shQuote(normalizePath(opt$outdir)), opt$image, basename(fp), basename(fn),
                   paste("--m", DBS, collapse = " "), tag)
    if (system(cmd) != 0) stop("SEA failed: ", tag)
    r <- read.delim(file.path(opt$outdir, paste0("sea_", tag), "sea.tsv"), comment.char = "#",
                    stringsAsFactors = FALSE)
    if (nrow(r)) runs[[tag]] <- r %>% mutate(set = nm, stage = st, direction = sub(".*_", "", nm),
                                           window = up, n_genes = length(sets[[nm]]))
  }
}
all <- bind_rows(runs)
write.table(all, file.path(opt$outdir, "sea_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", nrow(all), " SEA rows; significant (q<0.05): ", sum(all$QVALUE < 0.05))
