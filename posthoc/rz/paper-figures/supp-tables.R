#!/usr/bin/env Rscript
# Supplementary gene tables, one workbook each, with a README sheet first:
#   Supplementary_File_1_DE_hlh-30_vs_WT.xlsx  DESeq2 hlh-30(tm1978) vs WT at each stage
#   Supplementary_File_2_kmeans_clusters.xlsx  k-means cluster of every clustered gene (Fig 3A, 3E, S6C)
#   Supplementary_File_3_compare_DE_lists.xlsx compare-DE gene lists that went into GO (Fig 4, S2, 9)
#   Supplementary_File_5_hlh-30_gene_evidence.xlsx  published studies behind the gene counts of Fig S7A and S8A
#   Supplementary_File_6_SEA_motifs.xlsx           every SEA motif with q < 0.05 (Fig S10)
#
# Stage tables: L4, Day 3 and Day 6 come from the 28-sample L4 to Day 6 set
# (the Fig 1 DEG counts 29/84/251), L1 from the 34-sample set (Fig 9). The same
# L4 to Day 6 contrasts from the 34-sample set (35/86/258 DEGs, used by Fig S5C
# and S7 to S11) are extra sheets.
# Cluster numbers are the paper numbers (kmeans-umap/ cache, renumbered by
# _common/kmeans-renumber.tsv). Compare-DE lists are compare-de/ output
# (--tested_in_both TRUE), the lists the GO panels use.
#
#   make -C compare-de RUN=VZRJXQ-ruv-rbe compare; make -C compare-de RUN=hlh30-all-ruv compare
#   aws s3 cp s3://celegans-pipeline-rayandrew/refs/gene_db.sqlite ../../../gene_db.sqlite
#   c_elegans.PRJNA13758.WS295.functional_descriptions.txt.gz from the WormBase WS295
#   release (species/c_elegans/PRJNA13758/annotation/) into ../../../.cache/wormbase/
#   Rscript supp-tables.R

suppressPackageStartupMessages({ library(optparse); library(openxlsx) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--compare", type = "character", default = "../compare-de"),
  make_option("--kmeans",  type = "character", default = "../kmeans-umap/cache"),
  make_option("--opposite", type = "character", default = "out/07-opposite-genes.csv"),
  make_option("--pub",     type = "character", default = "../hlh30-integration/published"),
  make_option("--sea",     type = "character", default = "../hlh30-integration/out/sea"),
  make_option("--gene_db", type = "character", default = "../../../gene_db.sqlite"),
  make_option("--descriptions", type = "character",
              default = "../../../.cache/wormbase/c_elegans.PRJNA13758.WS295.functional_descriptions.txt.gz"),
  make_option("--outdir",  type = "character", default = "out/supp")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

PADJ <- 0.05; LFC <- 1
# The pipeline DE tables give the gene ID as the name of about 9,600 genes
# that have only a sequence name, such as Y102A5C.6. Those names come from the
# WormBase annotation table of gene_db. A real pipeline name is kept, since
# gene_db names the 12 mitochondrial genes by sequence (MTCE.26 for ctc-1).
ann <- local({
  con <- DBI::dbConnect(RSQLite::SQLite(), opt$gene_db); on.exit(DBI::dbDisconnect(con))
  a <- DBI::dbGetQuery(con, "SELECT wormbase_gene, gene_name, sequence_name FROM annotations")
  setNames(ifelse(is.na(a$gene_name) | a$gene_name == "", a$sequence_name, a$gene_name), a$wormbase_gene)
})
wb_name <- function(id, name = id) {
  n <- unname(ann[id])
  ifelse(!is.na(name) & name != "" & !startsWith(name, "WBGene"), name, ifelse(is.na(n), id, n))
}
de_dir <- function(run) file.path(opt$compare, "cache", run, "de")
read_de <- function(run, file) read.csv(file.path(de_dir(run), file), stringsAsFactors = FALSE)

# The README sheet is a list of blocks, written top to bottom with a blank row
# between them: the bold title first, then text (one sentence per row), bold
# headings and small tables (one item per row, bold header).
title   <- function(x) list(kind = "title", x = x)
text    <- function(...) list(kind = "text", x = c(...))
heading <- function(x) list(kind = "heading", x = x)
tab     <- function(...) list(kind = "table", x = data.frame(..., check.names = FALSE))
cols    <- function(...) { v <- c(...); tab(Column = names(v), Description = unname(v)) }

write_readme <- function(wb, blocks) {
  addWorksheet(wb, "README")
  bold <- createStyle(textDecoration = "bold")
  r <- 1
  for (b in blocks) {
    if (b$kind == "table") {
      writeData(wb, "README", b$x, startRow = r, headerStyle = bold)
      r <- r + nrow(b$x) + 2
    } else {
      writeData(wb, "README", b$x, startRow = r, colNames = FALSE)
      if (b$kind != "text") addStyle(wb, "README", bold, rows = r, cols = 1)
      # a heading sits right on top of its table
      r <- r + length(b$x) + (b$kind != "heading")
    }
  }
  setColWidths(wb, "README", 1:4, c(32, 70, 30, 30))
}

# openxlsx has no rich text, so each README string that names hlh-30 (with its
# allele) is rewritten in the shared strings as runs with the name in italics.
# The title is bold, so all its runs are bold too. Strings are shared between
# sheets, so a data cell that holds the same text is pointed at a plain copy.
GENE <- "hlh-30(\\(tm1978\\))?"
italicize <- function(wb, n_readme, title) {
  done <- integer(0)
  for (i in seq_len(n_readme)) {
    si <- wb$sharedStrings[[i]]
    if (!grepl(GENE, si)) next
    txt <- sub('^<si><t xml:space="preserve">(.*)</t></si>$', "\\1", si)
    b <- if (identical(txt, gsub("&", "&amp;", title, fixed = TRUE))) "<b/>" else ""
    run <- function(x, it) if (nzchar(x)) sprintf('<r>%s<t xml:space="preserve">%s</t></r>',
      if (it || nzchar(b)) sprintf("<rPr>%s%s</rPr>", b, if (it) "<i/>" else "") else "", x) else ""
    m <- gregexpr(GENE, txt)[[1]]
    plain <- substring(txt, c(1, m + attr(m, "match.length")), c(m - 1, nchar(txt)))
    it <- regmatches(txt, list(m))[[1]]
    body <- c(rbind(vapply(plain[-length(plain)], run, "", it = FALSE), vapply(it, run, "", it = TRUE)),
              run(plain[length(plain)], FALSE))
    wb$sharedStrings[[i]] <- paste0("<si>", paste(body, collapse = ""), "</si>")
    wb$sharedStrings[[length(wb$sharedStrings) + 1]] <- si
    done[as.character(i - 1)] <- length(wb$sharedStrings) - 1
  }
  for (k in seq_along(wb$worksheets)[-1]) {
    sd <- wb$worksheets[[k]]$sheet_data
    hit <- which(sd$t == 1 & sd$v %in% names(done))
    if (length(hit)) sd$v[hit] <- as.character(done[sd$v[hit]])
  }
  attr(wb$sharedStrings, "uniqueCount") <- length(wb$sharedStrings)
}

save_book <- function(sheets, readme, file) {
  wb <- createWorkbook()
  write_readme(wb, readme)
  n_readme <- length(wb$sharedStrings)
  bold <- createStyle(textDecoration = "bold")
  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets[[nm]], headerStyle = bold)
    freezePane(wb, nm, firstRow = TRUE)
    setColWidths(wb, nm, seq_along(sheets[[nm]]), "auto")
  }
  italicize(wb, n_readme, readme[[1]]$x)
  saveWorkbook(wb, file.path(opt$outdir, file), overwrite = TRUE)
  message(file, ": ", paste(sprintf("%s (%d)", names(sheets), sapply(sheets, nrow)), collapse = ", "))
}

# ---- 1. DESeq2 hlh-30 vs WT per stage ----------------------------------------
de_sheet <- function(run, file) {
  d <- read_de(run, file)
  d$direction <- ifelse(is.na(d$padj), "filtered",
                 ifelse(d$padj < PADJ & d$log2FoldChange >=  LFC, "up",
                 ifelse(d$padj < PADJ & d$log2FoldChange <= -LFC, "down", "not DEG")))
  d <- d[order(d$direction == "filtered", d$padj, d$pvalue), ]
  data.frame(gene_id = d$gene_id, gene_name = wb_name(d$gene_id, d$gene_name), baseMean = d$baseMean,
             log2FC = d$log2FoldChange, lfcSE = d$lfcSE, pvalue = d$pvalue, padj = d$padj,
             direction = d$direction, check.names = FALSE)
}
de <- list(
  "L1"               = de_sheet("hlh30-all-ruv",  "tm1978_L1_vs_WT_L1.csv"),
  "L4"               = de_sheet("VZRJXQ-ruv-rbe", "hlh-30_L4_vs_N2_L4.csv"),
  "Day 3"            = de_sheet("VZRJXQ-ruv-rbe", "hlh-30_Day3_vs_N2_Day3.csv"),
  "Day 6"            = de_sheet("VZRJXQ-ruv-rbe", "hlh-30_Day6_vs_N2_Day6.csv"),
  "L4 (34-sample set)"    = de_sheet("hlh30-all-ruv",  "tm1978_L4_vs_WT_L4.csv"),
  "Day 3 (34-sample set)" = de_sheet("hlh30-all-ruv",  "tm1978_Day3_vs_WT_Day3.csv"),
  "Day 6 (34-sample set)" = de_sheet("hlh30-all-ruv",  "tm1978_Day6_vs_WT_Day6.csv"))
n_of <- function(d, x) sum(d$direction == x)
set_of <- c("34-sample set", rep("28-sample set", 3), rep("34-sample set", 3))
used_in <- c("Fig. 9B to 9D, S5C and S7 to S11", rep("Fig. 1E, 1F and S1A", 3), rep("Fig. S5C and S7 to S11", 3))
save_book(de, list(
  title("Supplementary File 1. Differential expression of hlh-30(tm1978) vs WT at L1, L4, Day 3 and Day 6"),
  text("Each sheet lists every gene tested in one contrast of hlh-30(tm1978) with wild type (WT) at one stage.",
       sprintf("A gene is a differentially expressed gene (DEG) when padj < %g and |log2FC| >= %g.", PADJ, LFC),
       "Column 'direction' marks the DEGs, where:",
       "Up means higher in hlh-30(tm1978) than in WT, and down means lower in hlh-30(tm1978) than in WT.",
       "Two DESeq2 models were fit, one on the 28-sample set (L4, Day 3 and Day 6) and one on the 34-sample set (L1, L4, Day 3 and Day 6), and the Sheets table gives the model of each sheet.",
       "L4, Day 3 and Day 6 are in both models, so each of these stages has two sheets with slightly different DEG counts.",
       "The k-means clusters (Supplementary File 2) and the age comparisons (Supplementary File 3) use other contrasts and are not in this file."),
  heading("Sheets"),
  tab(Sheet = names(de), `Sample set` = set_of, `Used in` = used_in,
      `DEGs (up, down)` = sapply(de, function(d) sprintf("%d (%d, %d)", n_of(d, "up") + n_of(d, "down"),
                                                       n_of(d, "up"), n_of(d, "down")))),
  heading("Columns"),
  cols(gene_id   = "WormBase gene ID",
       gene_name = "WormBase gene name (the sequence name when the gene has no other name)",
       baseMean  = "Mean normalized count over all samples of the set",
       log2FC    = "log2 fold change of hlh-30(tm1978) over WT, shrunk with ashr (positive means higher in the mutant)",
       lfcSE     = "Posterior standard deviation of the shrunken log2FC",
       pvalue    = "Wald test P value (from the log2FC before shrinkage)",
       padj      = "Benjamini-Hochberg adjusted P value (empty when DESeq2 independent filtering removed the gene because of its low mean count, so it cannot be a DEG)",
       direction = "up or down for a DEG, not DEG for any other gene with a padj, filtered when padj is empty")),
  "Supplementary_File_1_DE_hlh-30_vs_WT.xlsx")

# ---- 2. k-means clusters ------------------------------------------------------
km_sheet <- function(run) {
  k <- read.csv(file.path(opt$kmeans, run, "cluster_assignments.csv"), stringsAsFactors = FALSE)
  k$gene_name <- wb_name(k$gene_id, k$gene_name)
  k <- k[order(k$cluster, k$gene_name), ]
  data.frame(gene_id = k$gene_id, gene_name = k$gene_name, cluster = k$cluster)
}
km <- list("Fig. 3A (L4 to Day 6)" = km_sheet("VZRJXQ-ruv-rbe"),
           "Fig. 3E (Day 3 and Day 6)" = km_sheet("VZRJXQ_no_l4-ruv-rbe"),
           "Fig. S6C (all stages)" = km_sheet("hlh30-all-ruv"))
# k-means clustered the union of the DEG lists of each run, one list per
# comparison. One row per comparison, one column per sheet.
STAGES <- c("L1", "L4", "Day3", "Day6")
de_lists <- function(run) {
  f <- list.files(file.path(opt$kmeans, run, "de-genes"), "_genes_de\\.txt$", full.names = TRUE)
  setNames(lapply(f, function(x) unique(readLines(x))), sub("_genes_de\\.txt$", "", basename(f)))
}
comparison <- function(id) {
  p <- strsplit(strsplit(id, "_vs_")[[1]], "_")
  g <- sub("^(N2|WT)$", "WT", sub("^(hlh-30|tm1978)$", "hlh-30", sapply(p, `[`, 1)))
  s <- sapply(p, `[`, 2); stage <- function(x) sub("Day", "Day ", x)
  if (g[1] != g[2]) return(sprintf("hlh-30 vs WT at %s", stage(s[1])))
  s <- s[order(match(s, STAGES), decreasing = TRUE)]
  sprintf("%s, %s vs %s", g[1], stage(s[1]), stage(s[2]))
}
runs <- c("Fig. 3A (L4 to Day 6)" = "VZRJXQ-ruv-rbe", "Fig. 3E (Day 3 and Day 6)" = "VZRJXQ_no_l4-ruv-rbe",
          "Fig. S6C (all stages)" = "hlh30-all-ruv")
degs <- lapply(runs, de_lists)
for (s in names(runs)) stopifnot(setequal(unique(unlist(degs[[s]])), km[[s]]$gene_id))
n_deg <- lapply(degs, function(l) setNames(lengths(l), sapply(names(l), comparison)))
cmp <- unique(unlist(lapply(n_deg, names)))
# rows in stage order: hlh-30 vs WT, then the stage pairs of WT, then of hlh-30
st <- sub("Day", "Day ", STAGES)
pairs <- unlist(lapply(2:4, function(i) sprintf("%s vs %s", st[i], st[seq_len(i - 1)])))
cmp <- intersect(c(sprintf("hlh-30 vs WT at %s", st), paste0("WT, ", pairs), paste0("hlh-30, ", pairs)), cmp)
stopifnot(length(cmp) == length(unique(unlist(lapply(n_deg, names)))))
deg_tab <- data.frame(Contrast = cmp, lapply(n_deg, function(n) unname(n[cmp])), check.names = FALSE)

# one row per sheet, one column per cluster
sizes <- do.call(rbind, lapply(km, function(d) {
  n <- table(factor(d$cluster, levels = 1:6))
  setNames(as.list(ifelse(n == 0, NA, n)), paste("Cluster", 1:6))
}))
save_book(km, list(
  title("Supplementary File 2. Genes of each k-means cluster"),
  text("Each sheet lists the DEGs that were clustered for one figure panel, with the k-means cluster of each gene.",
       "Cluster numbers are those shown in the figures.",
       "k-means grouped the genes by their mean VST expression in each genotype and stage, z-scored per gene.",
       "Each sample set has its own DESeq2 model, and a gene is a DEG when padj < 0.05 and |log2FC| >= 1 in at least one contrast of that model.",
       "The contrasts are hlh-30(tm1978) vs WT at each stage and each pair of stages within each genotype.",
       "A gene that is a DEG in more than one contrast is in the sheet only once, so the number of genes in a sheet is smaller than the sum of its DEG counts.",
       "The Fig. 3E model has only the Day 3 and Day 6 samples, so its DEG counts differ from those of the Fig. 3A model and of Supplementary File 1."),
  heading("Sheets"),
  tab(Sheet = names(km),
      `Sample set` = c("28-sample set (L4, Day 3 and Day 6)", "22-sample set (Day 3 and Day 6 samples of the 28-sample set)",
                       "34-sample set (L1, L4, Day 3 and Day 6)"),
      `k (clusters)` = c(6, 4, 6), Genes = sapply(km, nrow)),
  heading("DEGs per contrast"),
  tab(deg_tab),
  heading("Genes per cluster"),
  tab(Sheet = names(km), as.data.frame(lapply(as.data.frame(sizes), unlist), check.names = FALSE)),
  heading("Columns"),
  cols(gene_id = "WormBase gene ID", gene_name = "WormBase gene name (the sequence name when the gene has no other name)", cluster = "k-means cluster number")),
  "Supplementary_File_2_kmeans_clusters.xlsx")

# ---- 3. compare-DE gene lists -------------------------------------------------
# a = mutant, b = WT in every pair; lfc_diff = a - b
LISTS <- c("WT exclusive", "hlh-30 exclusive", "differential", "shared")
DIRS  <- c("up", "down", "more up in WT", "more up in hlh-30", "more down in WT", "more down in hlh-30",
           "up in WT, down in hlh-30", "down in WT, up in hlh-30")
# names as in Supplementary File 1, which covers every gene tested in both runs
de_names <- unlist(lapply(unname(de), function(d) setNames(d$gene_name, d$gene_id)))
de_names <- de_names[!duplicated(names(de_names))]
cmp_sheet <- function(run, pair, mut, wt) {
  m <- read.csv(file.path(opt$compare, "out", run, paste0(pair, "_compared_genes.csv")),
                stringsAsFactors = FALSE, check.names = FALSE)
  m <- m[m$category != "not_de", ]
  lm <- m[[paste0("log2FoldChange.", mut)]]; lw <- m[[paste0("log2FoldChange.", wt)]]
  dir <- function(x) ifelse(x > 0, "up", "down")
  m$list <- ifelse(m$category == paste0(wt, "_exclusive"), "WT exclusive",
            ifelse(m$category == paste0(mut, "_exclusive"), "hlh-30 exclusive",
            ifelse(m$category == "shared", "shared", "differential")))
  m$direction <- ifelse(m$list == "WT exclusive", dir(lw),
                 ifelse(m$list == "hlh-30 exclusive", dir(lm),
                 ifelse(m$list == "shared", dir(lw),
                 # differential: by the signs of the two age changes, then by size
                 ifelse(lw > 0 & lm > 0, ifelse(lw > lm, "more up in WT", "more up in hlh-30"),
                 ifelse(lw < 0 & lm < 0, ifelse(lw < lm, "more down in WT", "more down in hlh-30"),
                 ifelse(lw > 0, "up in WT, down in hlh-30", "down in WT, up in hlh-30"))))))
  m$gene_name <- wb_name(m$gene_id, unname(de_names[m$gene_id]))
  m <- m[order(match(m$list, LISTS), match(m$direction, DIRS), m$gene_name), ]
  data.frame(gene_id = m$gene_id, gene_name = m$gene_name, gene_list = m$list, direction = m$direction,
             log2FC_WT = m[[paste0("log2FoldChange.", wt)]], padj_WT = m[[paste0("padj.", wt)]],
             log2FC_hlh30 = m[[paste0("log2FoldChange.", mut)]], padj_hlh30 = m[[paste0("padj.", mut)]],
             log2FC_difference = -m$lfc_diff, check.names = FALSE)
}
V <- "VZRJXQ-ruv-rbe"
cmp <- list(
  "L1 to L4"    = cmp_sheet("hlh30-all-ruv", "tm1978_L4_vs_tm1978_L1__X__WT_L4_vs_WT_L1",
                            "tm1978_L4_vs_tm1978_L1", "WT_L4_vs_WT_L1"),
  "L4 to Day 3" = cmp_sheet(V, "hlh-30_Day3_vs_hlh-30_L4__X__N2_Day3_vs_N2_L4",
                            "hlh-30_Day3_vs_hlh-30_L4", "N2_Day3_vs_N2_L4"),
  "L4 to Day 6" = cmp_sheet(V, "hlh-30_Day6_vs_hlh-30_L4__X__N2_Day6_vs_N2_L4",
                            "hlh-30_Day6_vs_hlh-30_L4", "N2_Day6_vs_N2_L4"),
  "Day 3 to Day 6" = cmp_sheet(V, "hlh-30_Day6_vs_hlh-30_Day3__X__N2_Day6_vs_N2_Day3",
                               "hlh-30_Day6_vs_hlh-30_Day3", "N2_Day6_vs_N2_Day3"),
  # the existing new-analyses list (fc-quadrant/), kept as it is
  "Opposite direction" = read.csv(opt$opposite, stringsAsFactors = FALSE, check.names = FALSE))
opp <- cmp[["Opposite direction"]]
cmp <- cmp[names(cmp) != "Opposite direction"]
# the fc-quadrant list writes D3/D6 and rounds; take full values from the transition sheets
opp$transition <- sub("\\bD([36])\\b", "Day \\1", opp$transition)
stopifnot(opp$transition %in% names(cmp))
opp <- do.call(rbind, lapply(seq_len(nrow(opp)), function(i) {
  m <- cmp[[opp$transition[i]]]; m <- m[m$gene_id == opp$gene_id[i], ]
  stopifnot(nrow(m) == 1, m$gene_list == "differential", sign(m$log2FC_WT) != sign(m$log2FC_hlh30))
  data.frame(transition = opp$transition[i], m[, c("gene_id", "gene_name", "log2FC_WT", "padj_WT", "log2FC_hlh30", "padj_hlh30")],
             check.names = FALSE)
}))
# Most genes have no curated WormBase description (the pipeline tables carry
# only that one), so fall back to the WormBase automated description.
desc <- local({
  x <- readLines(gzfile(opt$descriptions))
  x <- x[!startsWith(x, "#")]
  rec <- split(x, cumsum(x == "="))
  field <- function(r, key) {
    i <- grep(paste0("^", key, ":"), r)
    if (!length(i)) return(NA_character_)
    j <- i
    while (j < length(r) && !grepl("^[A-Z][a-z]+( [a-z]+)? description:|^=$", r[j + 1])) j <- j + 1
    v <- trimws(gsub("\\s+", " ", sub(paste0("^", key, ":"), "", paste(r[i:j], collapse = " "))))
    if (v %in% c("none available", "not known", "")) NA_character_ else v
  }
  rec <- lapply(rec, function(r) r[r != "="])
  rec <- rec[vapply(rec, function(r) length(r) > 0 && startsWith(r[1], "WBGene"), logical(1))]
  data.frame(gene_id = sub("\t.*", "", vapply(rec, `[`, "", 1)),
             concise = vapply(rec, field, "", "Concise description"),
             automated = vapply(rec, field, "", "Automated description"))
})
d <- desc[match(opp$gene_id, desc$gene_id), ]
stopifnot(!is.na(d$gene_id))
opp$description <- ifelse(!is.na(d$concise), d$concise, d$automated)
opp$description_source <- ifelse(!is.na(d$concise), "curated", ifelse(!is.na(d$automated), "automated", NA))
# one row per list and direction, one column per transition
combos <- unique(do.call(rbind, lapply(cmp, function(d) d[, c("gene_list", "direction")])))
combos <- combos[order(match(combos$gene_list, LISTS), match(combos$direction, DIRS)), ]
counts <- data.frame(`Gene list` = combos$gene_list, Direction = combos$direction, check.names = FALSE)
for (s in names(cmp)) counts[[s]] <- sapply(seq_len(nrow(combos)), function(i)
  sum(cmp[[s]]$gene_list == combos$gene_list[i] & cmp[[s]]$direction == combos$direction[i]))
save_book(c(cmp, list("Opposite direction" = opp)), list(
  title("Supplementary File 3. Comparative differential expression gene lists, with the opposite-direction genes"),
  text("Each transition sheet lists the genes that change with age in WT, in hlh-30(tm1978) or in both.",
       "For each transition, the age change in each genotype is one contrast of the DESeq2 model of the sample set (see Sheets).",
       sprintf("A gene is a DEG in a genotype when padj < %g and |log2FC| >= %g.", PADJ, LFC),
       "Only genes with a padj in both genotypes are compared.",
       "Fig. 4, S2 and 9F to 9I show the GO enrichment of selected gene lists."),
  heading("Gene lists"),
  tab(`Gene list` = c("WT exclusive", "hlh-30 exclusive", "differential", "shared"),
      Direction = c("up or down", "up or down", "see Directions of the differential genes", "up or down"),
      Meaning = c("DEG in WT only, up or down with age",
                  "DEG in hlh-30(tm1978) only, up or down with age",
                  "DEG in both genotypes, and the two log2FC values differ by 1 or more",
                  "DEG in both genotypes, and the two log2FC values differ by less than 1, so both change in the same direction")),
  heading("Directions of the differential genes"),
  tab(Direction = DIRS[3:8],
      Meaning = c("Both genotypes go up, and WT goes up more", "Both genotypes go up, and hlh-30(tm1978) goes up more",
                  "Both genotypes go down, and WT goes down more", "Both genotypes go down, and hlh-30(tm1978) goes down more",
                  "WT goes up, and hlh-30(tm1978) goes down", "WT goes down, and hlh-30(tm1978) goes up")),
  text("Fig. 9F and 9G show the GO enrichment of the L1 to L4 genes labeled 'more up in WT' or 'up in WT, down in hlh-30'."),
  heading("Sheets"),
  tab(Sheet = c(names(cmp), "Opposite direction"),
      `Sample set` = c("34-sample set", rep("28-sample set", 3), "As in its transition sheet"),
      Genes = c(sapply(cmp, nrow), nrow(opp))),
  text(sprintf("The Opposite direction sheet lists the %d genes that are DEGs in both genotypes but change in opposite directions with age.", nrow(opp)),
       "They are the orange points of the quadrant plots in Fig. 4D and 9E, and they are also in the differential lists of the transition sheets."),
  heading("Genes per list"),
  tab(counts),
  heading("Columns of the transition sheets"),
  cols(gene_id = "WormBase gene ID", gene_name = "WormBase gene name (the sequence name when the gene has no other name)",
       gene_list = "WT exclusive, hlh-30 exclusive, differential or shared (see Gene lists)",
       direction = "up (higher at the later stage) or down (lower at the later stage), or one of the Directions of the differential genes",
       log2FC_WT = "log2 fold change of the later stage over the earlier stage in WT",
       padj_WT = "Adjusted P value of the age change in WT",
       log2FC_hlh30 = "log2 fold change of the later stage over the earlier stage in hlh-30(tm1978)",
       padj_hlh30 = "Adjusted P value of the age change in hlh-30(tm1978)",
       log2FC_difference = "log2FC_WT minus log2FC_hlh30"),
  heading("Columns of the Opposite direction sheet"),
  cols(transition = "The age transition (the name of its transition sheet)", gene_id = "WormBase gene ID", gene_name = "WormBase gene name (the sequence name when the gene has no other name)",
       "log2FC_WT, padj_WT" = "Age change in WT", "log2FC_hlh30, padj_hlh30" = "Age change in hlh-30(tm1978)",
       description = "WormBase WS295 gene description: the curated description when the gene has one, otherwise the automated description",
       description_source = "curated (written by WormBase curators) or automated (generated by WormBase from GO, expression and orthology data)")),
  "Supplementary_File_3_compare_DE_lists.xlsx")

# ---- 5. published evidence for the hlh-30-linked autophagy and lysosome genes --
# The number after each gene name in Fig S7A and S8A is n_papers_expr_reduced_hlh30
# (published-overlap.R); the other evidence types are listed but not counted.
ev  <- read.delim(file.path(opt$pub, "new_autophagy", "hlh30_autophagy_evidence.tsv"),
                  stringsAsFactors = FALSE, quote = "", na.strings = "")
gs  <- read.delim(file.path(opt$pub, "new_autophagy", "hlh30_autophagy_genes_summary.tsv"),
                  stringsAsFactors = FALSE, quote = "", na.strings = "")
lys <- read.delim(file.path(opt$pub, "new_lysosome", "sun2020_lysosome_genes.tsv"),
                  stringsAsFactors = FALSE, quote = "")
EVID <- c(expr_reduced_hlh30 = "lower without hlh-30", induced_hlh30_oe = "higher with hlh-30 overexpression",
          chip_binding = "HLH-30 ChIP binding", motif_only = "promoter motif only",
          not_dependent = "not dependent on hlh-30")
CLASS <- c(core_machinery = "core machinery", maturation_fusion = "maturation / fusion",
           selective_autophagy = "selective autophagy", lysosome_or_other = "lysosome or other")
is_true <- function(x) x %in% c("True", "TRUE", TRUE)
counted <- gs$n_papers_expr_reduced_hlh30 > 0
shown <- trimws(paste(ifelse(counted & is_true(gs$in_core_list), "Fig. S7A", ""),
                      ifelse(counted & gs$gene_id %in% lys$gene_id, "Fig. S8A", "")))
genes <- data.frame(gene_id = gs$gene_id, gene_name = gs$gene_name, gene_class = CLASS[gs$gene_class],
                    shown_in = shown, n_studies_lower = gs$n_papers_expr_reduced_hlh30,
                    studies_lower = gs$papers_expr_reduced_hlh30,
                    studies_higher_with_overexpression = gs$papers_induced_hlh30_oe,
                    studies_chip_binding = gs$papers_chip_binding, studies_motif_only = gs$papers_motif_only,
                    studies_not_dependent = gs$papers_not_dependent, check.names = FALSE)
genes <- genes[order(genes$shown_in == "", -genes$n_studies_lower, genes$gene_name), ]
studies <- data.frame(gene_id = ev$gene_id, gene_name = ev$gene_name, name_in_study = ev$source_gene,
                      gene_class = CLASS[ev$gene_class], evidence = EVID[ev$evidence],
                      counted_in_figure = ifelse(ev$evidence == "expr_reduced_hlh30", "yes", "no"),
                      hlh30_perturbation = ev$hlh30_perturbation, condition = ev$condition,
                      stage_tissue = ev$stage_tissue, assay = ev$assay, result = ev$result,
                      study = ev$paper, figure = ev$figure, pmcid = ev$pmcid, doi = ev$doi,
                      quote = ev$quote, notes = ev$notes, check.names = FALSE)
studies <- studies[order(studies$gene_name, studies$study, na.last = TRUE), ]
save_book(list("Genes" = genes, "Studies" = studies), list(
  title("Supplementary File 5. Studies that link each autophagy and lysosome gene to hlh-30"),
  text("This file lists the published studies behind the number after each gene name in Fig. S7A and S8A (column n_studies_lower).",
       "A study counts for a gene when the gene's mRNA or transcriptional reporter was lower, or its induction was lost or reduced, in an hlh-30 mutant or after hlh-30 RNAi.",
       "Fig. S7A shows the genes of the core autophagy list (Zhang and Melendez 2025, Table 1) with at least one counted study.",
       "Fig. S8A shows the genes of the lysosome list (Sun et al. 2020, Supplementary file 1) with at least one counted study.",
       "The other evidence types are listed but not counted.",
       sprintf("The file has %d genes and %d studies, and %d genes have at least one counted study.",
               nrow(genes), length(unique(ev$paper)), sum(counted))),
  heading("Sheets"),
  tab(Sheet = c("Genes", "Studies"), Rows = c(nrow(genes), nrow(studies)),
      Description = c("One row per gene", "One row per result, so a study can have several rows for one gene (such as one per condition)")),
  heading("Evidence types"),
  tab(Evidence = unname(EVID), `Counted in the figures` = c("yes", "no", "no", "no", "no"),
      Meaning = c("The gene was lower, or its induction was lost or reduced, without hlh-30",
                  "The gene was higher when HLH-30 was overexpressed",
                  "The study showed HLH-30 ChIP binding at the promoter",
                  "The study predicted a CLEAR, E-box or HLH-30 site in the promoter, with no binding data",
                  "The study removed hlh-30 and the change in the gene was still there")),
  heading("Columns of the Genes sheet"),
  cols(gene_id = "WormBase gene ID", gene_name = "WormBase gene name",
       gene_class = "'core machinery', 'maturation / fusion', 'selective autophagy' or 'lysosome or other'",
       shown_in = "The figure panel that shows the gene (empty when the gene has no counted study or is not in the core autophagy or lysosome list)",
       n_studies_lower = "Number of studies in which the gene was lower without hlh-30",
       studies_lower = "The studies in which the gene was lower without hlh-30",
       studies_higher_with_overexpression = "The studies in which the gene was higher with hlh-30 overexpression",
       studies_chip_binding = "The studies that show HLH-30 binding at the gene",
       studies_motif_only = "The studies that predict a promoter motif only",
       studies_not_dependent = "The studies that found the gene does not depend on hlh-30"),
  heading("Columns of the Studies sheet"),
  cols(gene_id = "WormBase gene ID (empty when the study's gene name matches more than one gene, such as ctsa)",
       gene_name = "WormBase gene name", name_in_study = "The gene name as written in the study",
       gene_class = "As in the Genes sheet", evidence = "The evidence type (see Evidence types)",
       counted_in_figure = "yes for the evidence type that is counted in the figures",
       hlh30_perturbation = "How hlh-30 was changed, such as a mutant allele or RNAi",
       condition = "The condition of the experiment", stage_tissue = "The stage and tissue",
       assay = "The assay", result = "The result",
       study = "First author, year and journal", figure = "The figure or table of the study that shows the result",
       pmcid = "PubMed Central ID", doi = "DOI of the study",
       quote = "Supporting text from the study", notes = "Curation notes")),
  "Supplementary_File_5_hlh-30_gene_evidence.xlsx")

# ---- 6. SEA motif enrichment ----------------------------------------------------
sea <- read.delim(file.path(opt$sea, "sea_significant.tsv"), stringsAsFactors = FALSE, quote = "",
                  check.names = FALSE)
grid <- read.delim(file.path(opt$sea, "sea_grid_factors.tsv"), stringsAsFactors = FALSE, quote = "")
SET <- c(L1 = "L1", L4 = "L4", D3 = "Day 3", D6 = "Day 6", any = "all stages")
st <- sub(" .*", "", sea$set)
DB <- ifelse(grepl("JASPAR", sea$DB), "JASPAR 2022 CORE nematodes", "CIS-BP 2.00 C. elegans")
motifs <- data.frame(stage = SET[st], direction = sub(".* ", "", sea$set),
                     promoter_bp = sea$window, database = DB, motif_id = sea$ID,
                     factor = sub(" \\(inf\\.\\)$", "", sea$tf), inferred = ifelse(grepl("\\(inf\\.\\)$", sea$tf), "yes", "no"),
                     consensus = sea$CONSENSUS, TP = sea$TP, TP_percent = sea[["TP."]],
                     FP = sea$FP, FP_percent = sea[["FP."]], enrichment_ratio = sea$ENR_RATIO,
                     q_value = sea$QVALUE,
                     factor_expressed = ifelse(sea$expressed %in% c(TRUE, "TRUE"), "yes", "no"),
                     in_figure = ifelse(sea$tf %in% grid$tf & sea$expressed %in% c(TRUE, "TRUE"), "yes", "no"),
                     check.names = FALSE)
save_book(list("Enriched motifs" = motifs), list(
  title("Supplementary File 6. All enriched SEA motifs"),
  text("This file lists every transcription factor motif enriched (SEA q < 0.05) in the promoters of the hlh-30(tm1978) vs WT DEGs of the 34-sample set (Fig. S10).",
       "The control promoters of a stage were those of the genes with padj >= 0.05 at that stage, and for the all stages sets those with padj >= 0.05 at all four stages.",
       "The test was SEA of MEME Suite 5.5.7.",
       "Fig. S10 shows a motif only when the gene of its factor is expressed at that stage (column factor_expressed).",
       sprintf("The file has %d enriched motifs, and %d of them have an expressed factor.",
               nrow(motifs), sum(motifs$factor_expressed == "yes"))),
  heading("Sheets"),
  tab(Sheet = "Enriched motifs", Rows = nrow(motifs), Description = "One row per enriched motif and DEG set"),
  heading("Columns"),
  cols(stage = "The stage of the DEG set (L1, L4, Day 3 or Day 6), or all stages for the union of the DEGs of the four stages (no motif was enriched at L4)",
       direction = "up or down, the DEGs higher or lower in hlh-30(tm1978) than in WT",
       promoter_bp = "The promoter length, 500 bp or 2,000 bp upstream of each TSS",
       database = "JASPAR 2022 CORE nematodes (non-redundant) or CIS-BP 2.00 C. elegans",
       motif_id = "The motif ID in its database",
       factor = "The transcription factor of the motif",
       inferred = "yes when the CIS-BP motif was inferred from a related DNA-binding domain",
       consensus = "The consensus sequence of the motif",
       TP = "Number of DEG promoters with a motif site",
       TP_percent = "Percent of DEG promoters with a motif site",
       FP = "Number of control promoters with a motif site",
       FP_percent = "Percent of control promoters with a motif site",
       enrichment_ratio = "SEA enrichment ratio (close to TP_percent divided by FP_percent)",
       q_value = "SEA q value",
       factor_expressed = "yes when the factor's gene has a mean DESeq2-normalized count of at least 10 over the samples of that stage (any stage for the all stages sets)",
       in_figure = sprintf("yes for the rows of the %d factors in Fig. S10, when the factor is expressed (the figure shows the smallest q of each factor for each DEG set and promoter length)", nrow(grid)))),
  "Supplementary_File_6_SEA_motifs.xlsx")
