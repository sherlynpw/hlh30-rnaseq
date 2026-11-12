# Export the GO:0006914 (autophagy) gene set, including descendant terms (GOALL),
# from org.Ce.eg.db into raw/go0006914_orgCeegdb.tsv. build_lists.py reads that file.
#   Rscript export_go.R

suppressMessages(library(org.Ce.eg.db))

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))

go <- AnnotationDbi::select(
  org.Ce.eg.db,
  keys = "GO:0006914", keytype = "GOALL",
  columns = c("WORMBASE", "SYMBOL", "EVIDENCEALL")
)
go <- go[!is.na(go$WORMBASE), ]

collapse <- function(x) paste(sort(unique(x)), collapse = ",")
out <- do.call(rbind, lapply(split(go, go$WORMBASE), function(d) {
  data.frame(
    gene_id = d$WORMBASE[1],
    symbol = d$SYMBOL[1],
    evidence_codes = collapse(d$EVIDENCEALL),
    iea_only = all(d$EVIDENCEALL == "IEA")
  )
}))
out <- out[order(out$symbol), ]

info <- org.Ce.eg_dbInfo()
meta <- c(
  paste0("# org.Ce.eg.db ", packageVersion("org.Ce.eg.db")),
  paste0("# ", info$name, ": ", info$value)[info$name %in% c("GOSOURCEDATE", "EGSOURCEDATE")]
)
f <- file.path(here, "raw", "go0006914_orgCeegdb.tsv")
writeLines(meta, f)
suppressWarnings(write.table(out, f, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
cat(nrow(out), "genes written to", f, "\n")
