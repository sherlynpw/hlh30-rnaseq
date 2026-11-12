# Where each panel comes from

All paths are relative to `posthoc/rz/`. Every panel is one vector PDF plus a
PNG preview in `paper-figures/out/`. The panel letters follow the
Figures_Shorter_1 draft (Nov 2026). The staged files in `out/` and on S3 keep
the earlier labels in the "staged as" column (for example `fig-9j` is now Fig.
9F), and "not in the draft" marks a panel that the draft dropped. Pipeline
runs live under `s3://celegans-pipeline-rayandrew/results/<run>/`; run configs
are `samples/<run>.yaml` at the repo root.

Runs:

| run | design |
|---|---|
| `VZRJXQ-ruv-rbe` | WT + hlh-30(tm1978), L4 / Day 3 / Day 6, RUVs + removeBatchEffect |
| `VZRJXQ_no_l4-ruv-rbe` | same without L4 (Day 3 / Day 6 only) |
| `hlh30-all-ruv` | WT + hlh-30(tm1978), L1 / L4 / Day 3 / Day 6 |
| `pca-hlh30-all-ruv` | PCA-only run of the four-stage set |

## Figure 1

| panel | staged as | rendered file | script | run | pipeline input |
|---|---|---|---|---|---|
| 1B | fig-1b | `pca-stage/out/VZRJXQ-ruv-rbe/pca_stage.pdf` | `pca-stage/pca-stage.R` | VZRJXQ-ruv-rbe | `exploration/pca/` coordinates + variance |
| 1C | fig-1c | `sample-corr-heatmap/out/VZRJXQ-ruv-rbe/corr_corrected_pearson.pdf` | `sample-corr-heatmap/sample-corr-heatmap.R` | VZRJXQ-ruv-rbe | `deseq2/vst_counts.csv` (RUV-corrected), `counts/sample_meta.tsv` |
| 1D | fig-1d | `pca-genotype/out/VZRJXQ_no_l4-ruv-rbe/pca_genotype.pdf` | `pca-genotype/pca-genotype.R` | VZRJXQ_no_l4-ruv-rbe | `exploration/pca/` |
| 1E | fig-1e-day0 | `volcano/out/VZRJXQ-ruv-rbe/volcano_L4.pdf` | `volcano/render-volcano.R` | VZRJXQ-ruv-rbe | `deseq2/results/hlh-30_L4_vs_N2_L4.csv` |
| 1E | fig-1e-day3 | `volcano/out/VZRJXQ-ruv-rbe/volcano_Day3.pdf` | `volcano/render-volcano.R` | VZRJXQ-ruv-rbe | `deseq2/results/hlh-30_Day3_vs_N2_Day3.csv` |
| 1E | fig-1e-day6 | `volcano/out/VZRJXQ-ruv-rbe/volcano_Day6.pdf` | `volcano/render-volcano.R` | VZRJXQ-ruv-rbe | `deseq2/results/hlh-30_Day6_vs_N2_Day6.csv` |
| 1F | fig-1f | `gene-lists/out/fig-1f.pdf` | `gene-lists/gene-lists.R` | VZRJXQ-ruv-rbe | the three `hlh-30_<stage>_vs_N2_<stage>.csv`; DEGs at L4, D3 and D6 |

## Figure 3

| panel | staged as | rendered file | script | run | pipeline input |
|---|---|---|---|---|---|
| 3A | fig-3a | `kmeans-umap/out/fig-3a/kmeans_umap.pdf` | `kmeans-umap/kmeans-umap.R` | VZRJXQ-ruv-rbe | `deseq2/vst_counts.csv`, `clustering/kmeans/cluster_assignments.csv` (UMAP recomputed, seed 42) |
| 3B | fig-3b | `cluster-profiles/out/VZRJXQ-ruv-rbe/cluster_profiles.pdf` | `cluster-profiles/cluster-profiles.R` | VZRJXQ-ruv-rbe | `deseq2/vst_counts.csv`, `clustering/kmeans/cluster_assignments.csv`, `counts/sample_meta.tsv` |
| 3C | fig-3c | `go-plot/out-paper/VZRJXQ-ruv-rbe/[Fig 3C] GO-BP kmeans cluster 1.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `enrichment/go/kmeans/cluster_4_genes/go_BP.csv` (paper cluster 1) |
| 3D | fig-3d | `go-plot/out-paper/VZRJXQ-ruv-rbe/[Fig 3D] GO-CC kmeans cluster 1.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `enrichment/go/kmeans/cluster_4_genes/go_CC.csv` (paper cluster 1) |
| 3E | fig-3e | `kmeans-umap/out/fig-3e/kmeans_umap.pdf` | `kmeans-umap/kmeans-umap.R` | VZRJXQ_no_l4-ruv-rbe | as fig-3a |
| 3F | fig-3f | `cluster-profiles/out/VZRJXQ_no_l4-ruv-rbe/cluster_profiles.pdf` | `cluster-profiles/cluster-profiles.R` | VZRJXQ_no_l4-ruv-rbe | as fig-3b |
| 3G | fig-3g | `go-plot/out-paper/VZRJXQ_no_l4-ruv-rbe/[Fig 3G] GO-BP kmeans cluster 2.pdf` | `go-plot/render-dots.R` | VZRJXQ_no_l4-ruv-rbe | `enrichment/go/kmeans/cluster_1_genes/go_BP.csv` (paper cluster 2) |
| 3H | fig-3h | `go-plot/out-paper/VZRJXQ_no_l4-ruv-rbe/[Fig 3H] GO-BP kmeans cluster 3.pdf` | `go-plot/render-dots.R` | VZRJXQ_no_l4-ruv-rbe | `enrichment/go/kmeans/cluster_3_genes/go_BP.csv` (paper cluster 3) |

## Figure 4

Venns compare the WT and hlh-30 DEG lists of one transition, on the genes
DESeq2 tested in both genotypes (a padj in both contrasts). GO panels come
from compare-DE gene sets made the same way by `compare-de/` (`make
RUN=VZRJXQ-ruv-rbe`, `compare_de.R --tested_in_both TRUE`), with GO under
`enrichment/go/compare_de_tested_both/<set>/`. The pipeline's own
`enrichment/go/compare_de/` (no common filter) is not used.

| panel | staged as | rendered file | script | run | pipeline input |
|---|---|---|---|---|---|
| not in the draft | fig-4a | `venn-degs/out/shared/venn_B_l4_d3.pdf` | `venn-degs/venn-degs.R` | VZRJXQ-ruv-rbe | `deseq2/results/N2_Day3_vs_N2_L4.csv`, `hlh-30_L4_vs_hlh-30_Day3.csv` |
| 4A | fig-4a-options | `fc-quadrant/out/quadrant_B_l4_d3.pdf` | `fc-quadrant/fc-quadrant.R` | VZRJXQ-ruv-rbe | same two tables as 4A; x = log2FC WT, y = log2FC hlh-30, quadrant counts DEGs / all genes |
| 4C | fig-4b | `go-plot/out-paper/compare_DE/[Fig 4B] GO-CC L4-to-D3 exclusive-down WT.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `...N2_Day3_vs_N2_L4_exclusive_down/go_CC.csv` |
| 4B | fig-4c | `go-plot/out-paper/compare_DE/[Fig 4C] GO-CC L4-to-D3 exclusive-down hlh-30.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `...hlh-30_Day3_vs_hlh-30_L4_exclusive_down/go_CC.csv` |
| not in the draft | fig-4d | `venn-degs/out/shared/venn_C_l4_d6.pdf` | `venn-degs/venn-degs.R` | VZRJXQ-ruv-rbe | `N2_Day6_vs_N2_L4.csv`, `hlh-30_L4_vs_hlh-30_Day6.csv` |
| 4D | fig-4d-options | `fc-quadrant/out/quadrant_C_l4_d6.pdf` | `fc-quadrant/fc-quadrant.R` | VZRJXQ-ruv-rbe | same tables as 4D |
| 4E | fig-4e | `go-plot/out-paper/compare_DE/[Fig 4E] GO-BP L4-to-D6 exclusive-up hlh-30.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `...hlh-30_Day6_vs_hlh-30_L4_exclusive_up/go_BP.csv` |
| 4F | fig-4f | `go-plot/out-paper/compare_DE/[Fig 4F] GO-CC L4-to-D6 exclusive-down hlh-30.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `...hlh-30_Day6_vs_hlh-30_L4_exclusive_down/go_CC.csv` |
| not in the draft | fig-4g | `venn-degs/out/shared/venn_D_d3_d6.pdf` | `venn-degs/venn-degs.R` | VZRJXQ-ruv-rbe | `N2_Day6_vs_N2_Day3.csv`, `hlh-30_Day6_vs_hlh-30_Day3.csv` |
| 4G | fig-4g-options | `fc-quadrant/out/quadrant_D_d3_d6.pdf` | `fc-quadrant/fc-quadrant.R` | VZRJXQ-ruv-rbe | same tables as 4G |
| 4H | fig-4h | `go-plot/out-paper/compare_DE/[Fig 4H] GO-CC D3-to-D6 exclusive-down hlh-30.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `...hlh-30_Day6_vs_hlh-30_Day3_exclusive_down/go_CC.csv` |
| 4I | fig-4i | `go-plot/out-paper/compare_DE/[Fig 4I] GO-CC D3-to-D6 exclusive-up WT.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `...N2_Day6_vs_N2_Day3_exclusive_up/go_CC.csv` |

The full compare-DE set names are the `path` column of
`go-plot/targets-paper-figures.tsv`.

## Figure 8

| panel | staged as | rendered file | script | data |
|---|---|---|---|---|
| 8F | fig-8d | `chip-tracks/out/promoter_schematic.pdf` | `chip-tracks/promoter-schematic.R` | ce11 genome (BSgenome), GTF; CACGTG consensus scan of tts-1 and ztf-28 promoters |
| 8G | fig-8e | `chip-tracks/out/motif_counts.pdf` | `chip-tracks/motif-counts.R` | ce11 genome, GTF; consensus counts + JASPAR MA1449.1 for ~CAYRTG |

## Figure 9 and S6 (four-stage set)

| panel | staged as | rendered file | script | run | pipeline input |
|---|---|---|---|---|---|
| S6A | fig-9a | `pca-stage/out/pca-hlh30-all-ruv/pca_stage.pdf` | `pca-stage/pca-stage.R` | pca-hlh30-all-ruv | `exploration/pca/` |
| 9A | fig-9b | `pca-l1/out/pca_l1_genotype.pdf` | `pca-l1/pca-l1.R` | pca-hlh30-l1 | `exploration/pca/` (L1 samples only) |
| S6B | fig-9c | `sample-corr-heatmap/out/pca-hlh30-all-ruv/corr_corrected_pearson.pdf` | `sample-corr-heatmap/sample-corr-heatmap.R` | pca-hlh30-all-ruv | `deseq2/vst_counts.csv`, `counts/sample_meta.tsv` |
| 9B | fig-9d | `volcano/out/hlh30-all-ruv/volcano_L1.pdf` | `volcano/render-volcano.R` | hlh30-all-ruv | `deseq2/results/tm1978_L1_vs_WT_L1.csv` |
| 9C | fig-9e | `go-plot/out-paper/compare_DE/[Fig 9E] GO-BP L1 WT-vs-hlh-30 down.pdf` | `go-plot/render-dots.R` | hlh30-all-ruv | `enrichment/go/filter_de/filtered_tm1978_L1_vs_WT_L1_genes_down/go_BP.csv` |
| 9D | fig-9f | `go-plot/out-paper/compare_DE/[Fig 9F] GO-BP L1 WT-vs-hlh-30 up.pdf` | `go-plot/render-dots.R` | hlh30-all-ruv | `enrichment/go/filter_de/filtered_tm1978_L1_vs_WT_L1_genes_up/go_BP.csv` |
| S6C | fig-9g | `kmeans-umap/out/fig-9g/kmeans_umap.pdf` | `kmeans-umap/kmeans-umap.R` | hlh30-all-ruv | as fig-3a |
| S6D | fig-9h | `cluster-profiles/out/hlh30-all-ruv/cluster_profiles.pdf` | `cluster-profiles/cluster-profiles.R` | hlh30-all-ruv | as fig-3b |
| not in the draft | fig-9i | `venn-degs/out9/shared/venn_A_l1_l4.pdf` | `venn-degs/venn-degs.R` | hlh30-all-ruv | `WT_L1_vs_WT_L4.csv`, `tm1978_L4_vs_tm1978_L1.csv` |
| 9E | fig-9i-options | `fc-quadrant/out/quadrant_A_l1_l4.pdf` | `fc-quadrant/fc-quadrant.R` | hlh30-all-ruv | same tables as 9I |
| 9F | fig-9j | `go-plot/out-paper/compare_DE/[Fig 9J] GO-BP L1-to-L4 up-more-in-WT.pdf` | `go-plot/render-dots.R` | hlh30-all-ruv | `enrichment/go/compare_de/...genes_differential_up_more_in_WT_L4_vs_WT_L1_or_opposite/go_BP.csv` (genes up more in WT, or up in WT and down in hlh-30; `compare-de/split-differential.R`) |
| 9G | fig-9k | `go-plot/out-paper/compare_DE/[Fig 9K] GO-CC L1-to-L4 up-more-in-WT.pdf` | `go-plot/render-dots.R` | hlh30-all-ruv | `...genes_differential_up_more_in_WT_L4_vs_WT_L1_or_opposite/go_CC.csv` |
| 9H | fig-9l | `go-plot/out-paper/compare_DE/[Fig 9L] GO-BP L1-to-D0 exclusive-down hlh-30.pdf` | `go-plot/render-dots.R` | hlh30-all-ruv | `...genes_tm1978_L4_vs_tm1978_L1_exclusive_down/go_BP.csv` |
| 9I | fig-9m | `go-plot/out-paper/compare_DE/[Fig 9M] GO-CC L1-to-D0 exclusive-up WT.pdf` | `go-plot/render-dots.R` | hlh30-all-ruv | `...genes_WT_L4_vs_WT_L1_exclusive_up/go_CC.csv` |
| S6E | fig-s2f | `go-plot/out-paper/hlh30-all-ruv/[Fig S2F] GO-BP kmeans cluster 3.pdf` | `go-plot/render-dots.R` | hlh30-all-ruv | `enrichment/go/kmeans/cluster_3_genes/go_BP.csv` |

## Supplementary

| panel | staged as | rendered file | script | run | pipeline input |
|---|---|---|---|---|---|
| S1A | fig-s1a | `gene-lists/out/fig-s1a.pdf` | `gene-lists/gene-lists.R` | VZRJXQ-ruv-rbe | same three tables; DEGs at D3 and D6, not at L4 |
| S1B | fig-s1b | `kmeans-umap/out/fig-s1b/hdbscan_umap.pdf` | `kmeans-umap/kmeans-umap.R` | VZRJXQ-ruv-rbe | `clustering/hdbscan/cluster_assignments.csv`, `deseq2/vst_counts.csv` |
| S1C | fig-s1c | `pca-genotype/out/VZRJXQ_no_l4-ruv-rbe/pca_cohort.pdf` | `pca-genotype/pca-genotype.R --color-by cohort` | VZRJXQ_no_l4-ruv-rbe | the fig 1D coordinates, coloured by the sample sheet's `temperature` (20C or RT = 22 C); the cohort test behind it is `cohort-check/` |
| S2A | fig-s2a | `go-plot/out-paper/compare_DE/[Fig S2A] GO-BP L4-to-D3 exclusive-up WT.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `...N2_Day3_vs_N2_L4_exclusive_up/go_BP.csv` |
| S2B | fig-s2b | `go-plot/out-paper/compare_DE/[Fig S2B] GO-BP L4-to-D3 exclusive-down WT.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `...N2_Day3_vs_N2_L4_exclusive_down/go_BP.csv` |
| S2C | fig-s2c | `go-plot/out-paper/compare_DE/[Fig S2C] GO-BP L4-to-D6 exclusive-down hlh-30.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `...hlh-30_Day6_vs_hlh-30_L4_exclusive_down/go_BP.csv` |
| S2D | fig-s2d | `go-plot/out-paper/compare_DE/[Fig S2D] GO-BP L4-to-D6 exclusive-up WT.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `...N2_Day6_vs_N2_L4_exclusive_up/go_BP.csv` |
| S2E | fig-s2e | `go-plot/out-paper/compare_DE/[Fig S2E] GO-BP D3-to-D6 exclusive-down WT.pdf` | `go-plot/render-dots.R` | VZRJXQ-ruv-rbe | `...N2_Day6_vs_N2_Day3_exclusive_down/go_BP.csv` |
| S3A | fig-s3a | `germ-trajectory/out-s3a/germ_heatmap.pdf` | `germ-trajectory/germ-trajectory.R` | hlh30-all-ruv | `deseq2/vst_counts.csv`, `counts/sample_meta.tsv`; gene set built by `germ-trajectory/germ-genes.R` from VZRJXQ-ruv-rbe `enrichment/go/compare_de_tested_both/` (31 oocyte and germline genes) |
| S3B | fig-s3b | `germ-trajectory/out-s3a/sperm_heatmap.pdf` | `germ-trajectory/germ-trajectory.R` | hlh30-all-ruv | as fig-s3d |
| S3C | fig-s3c | `germ-trajectory/out-s3b/germ_meanline.pdf` | `germ-trajectory/germ-trajectory.R` | hlh30-all-ruv | as fig-s3a |
| S3D | fig-s3d | `germ-trajectory/out-s3c/sperm_meanline.pdf` | `germ-trajectory/germ-trajectory.R` | hlh30-all-ruv | as fig-s3a; sperm set = genes of the enriched pseudopodium term (18 genes) |
| S5A | fig-s5a | `chip-tracks/out/panelA_metagene.pdf` | `chip-tracks/chip-tracks.R` | ChIP-Atlas ce11 | HLH-30 bigWigs SRX331310/312/314/316 and matched inputs |
| S5B | fig-s5b | `chip-tracks/out/panelE_lin_style.pdf` | `chip-tracks/panel-e-lin.R` | ChIP-Atlas ce11 | same loci; one track per locus = mean of the 4 experiments (Lin et al. 2018 Fig 5e layout) |
| not in the draft | fig-s5b-lanes | `chip-tracks/out/panelE_tracks.pdf` | `chip-tracks/chip-tracks.R` | ChIP-Atlas ce11 | same bigWigs + MACS2 peak BEDs; loci ztf-28, tts-1, cpr-4 |
| S5C | fig-s5c | `motif-targets/out/hlh30-all-ruv/targets_dumbbell.pdf` | `motif-targets/render-targets.R` | hlh30-all-ruv | per-stage DEG lists + promoter scan for JASPAR MA1449.1 (canonical E-box only; BH over the 12 tests shown) |

## GO panels: hand-picked vs automatic terms

`go-plot/targets-paper-figures.tsv` records this in its `terms` column. A
non-empty cell lists the terms we chose by hand (with `=>` display renames);
an empty cell means the top terms by adjusted p, as the pipeline ranks them.

| hand-picked terms | automatic (top terms) |
|---|---|
| 3C, 3D (top 5 minus lytic vacuole, 2026-10-07), 3G, 3H, 4B, 4C, 4E, 4F, 4H, S2A to S2E, S6E, 9C, 9D, 9F, 9G, 9H, 9I | 4I |
