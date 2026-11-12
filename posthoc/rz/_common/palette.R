# One palette for every hlh-30 paper panel. Sourced, not run.
#
# One meaning per colour:
#   purple  = HLH-30 (the mutant genotype, its ChIP signal, motif, targets)
#   grey    = WT
#   red / blue = up / down (DEG direction; high / low in diverging scales)
#   stages and clusters = Okabe-Ito (colourblind-safe categorical set)
# Scales:
#   PAL$go_ramp   -log10 p in GO plots: blue (low) -> red (high)
#   PAL$hlh30_seq amounts of HLH-30 (sites, fold): white -> purple
#   PAL$div       diverging, RdBu: blue / white / red (z-score, correlation)

PAL <- list(
  hlh30 = "#7B3294",
  wt    = "#4D4D4D",
  up    = "#B2182B",
  down  = "#2166AC",
  ns    = "#BDBDBD",

  stage = c(L1 = "#009E73", L4 = "#E69F00", D3 = "#56B4E9", D6 = "#CC79A7"),

  # Okabe-Ito, in an order that keeps neighbouring clusters apart
  okabe = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7",
            "#0072B2", "#D55E00", "#F0E442", "#000000", "#999999"),

  # HLH-30 ChIP lanes: L4 in purples, late embryo in greys
  lane  = c("HLH-30  L4 rep1"  = "#7B3294", "HLH-30  L4 rep2"  = "#B892C9",
            "HLH-30  emb rep1" = "#4D4D4D", "HLH-30  emb rep2" = "#A6A6A6"),

  # blue -> red: Moreland's cool-warm map (bright, colourblind-safe)
  go_ramp   = c("#3B4CC0", "#6F92F3", "#A5C4F8", "#DDDDDD", "#F5A889", "#D85C4E", "#B40426"),
  # starts just off white, so a zero cell still shows as a tile
  hlh30_seq = c("#F3EDF5", "#DCC8E3", "#B892C9", "#7B3294", "#4C1D5C"),
  div       = c(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B")
)

# genotype colours keyed by every token the scripts use
PAL$geno <- c("WT" = PAL$wt, "N2" = PAL$wt,
              "hlh-30(tm1978)" = PAL$hlh30, "hlh-30" = PAL$hlh30, "tm1978" = PAL$hlh30)
