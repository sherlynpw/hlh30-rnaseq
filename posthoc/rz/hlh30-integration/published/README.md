# Published HLH-30 gene lists

Gene lists from published C. elegans HLH-30 studies, for overlap tests against our WT vs hlh-30(tm1978) DEGs.
Original supplementary files are in `raw/`.
`build_lists.py` turns them into the tidy TSVs in this folder:

```
uv run --with pandas --with openpyxl --with xlrd python build_lists.py
```

## TSV format

Every TSV has these columns:

- `gene_id`: WBGene ID, or empty when the name could not be mapped.
- `gene_name`: current WormBase name of `gene_id`, from `gene_db.sqlite` or the GTF.
- `source_gene`: the gene or sequence name exactly as the authors wrote it.
- `direction`: what membership means. `up` means induced in WT. `down` and `up` in the Grove/Harvald list refer to the hlh-30 mutant (see that section). `bound` means a ChIP peak at the promoter. `NA` means no direction (motif or score lists).
- `map_method`: how `gene_id` was found. `gene_name`, `sequence_name` or `other_names` means a case-insensitive match against `annotations` in `gene_db.sqlite`, tried in that order. `source` means the authors gave the WBGene ID. `gtf` means it came from the WBcel235.115 GTF.

Extra columns hold whatever scores the source gives.

"n mapped" below is the number of rows with a WBGene ID.
A few source names map to the same WBGene ID through `other_names`, so the number of distinct IDs can be smaller.

## Summary

| file | meaning | rows | mapped (distinct IDs) |
|---|---|---|---|
| visvikis2014_hlh30_dependent_saureus_induced.tsv | induced by S. aureus in WT, not in hlh-30 | 637 | 631 (631) |
| visvikis2014_hlh30_independent_saureus_induced.tsv | induced by S. aureus in WT and in hlh-30 | 188 | 184 (184) |
| visvikis2014_saureus_induced_wt.tsv | all S. aureus-induced genes in WT (S3 + S5) | 825 | 815 (815) |
| visvikis2014_hlh30_dependent_ebox.tsv | HLH-30-dependent induced genes with an E-box within 2 kb upstream | 70 | 70 (70) |
| orourke2013_hlh30_dependent_autophagy_qpcr.tsv | fasting-induced autophagy genes, not induced in hlh-30 (qPCR) | 6 | 6 (6) |
| lapierre2013_clear_motif_candidates.tsv | autophagy/lysosome genes with CLEAR-like sites (motif search only) | 17 | 17 (17) |
| harvald2017_known_hlh30_targets.tsv | previously reported hlh-30-regulated genes, apparently Grove 2009 microarray (122 down, 10 up in hlh-30) | 132 | 132 (128) |
| harvald2017_starvation_hlh30_protein_highB_clear.tsv | proteins whose starvation time course differs WT vs hlh-30 (B > 1.85) and whose gene has a CLEAR motif | 492 | 490 (490) |
| grove2009_hlh30_pbm_predicted_targets.tsv | genes predicted as HLH-30 targets from PBM motif plus conservation | 548 | 548 (548) |
| encode_modern_hlh30_l4_chip_promoter_bound.tsv | protein-coding genes with an HLH-30::GFP ChIP-seq summit at the promoter, L4 | 4017 | 4017 (4017) |

## Visvikis et al. 2014 (S. aureus infection)

Citation: Visvikis O, et al. Innate host defense requires TFEB-mediated transcription of cytoprotective and antimicrobial genes. Immunity 2014;40(6):896-909. doi:10.1016/j.immuni.2014.05.002. PMC4104614.

Source: `https://ars.els-cdn.com/content/image/1-s2.0-S1074761314001794-mmc2.xlsx` saved as `raw/1-s2.0-S1074761314001794-mmc2.xlsx`.
The Supplemental Information PDF (methods) is `raw/1-s2.0-S1074761314001794-mmc1.pdf`, from the same CDN as `-mmc1.pdf`.
The PMC copy of the author manuscript (NIHMS598830) sits behind a browser challenge, and Europe PMC answered "not open access".
The Elsevier CDN served the files without a browser.

Design: WT (N2) and hlh-30(tm1978) young adults on S. aureus or on non-pathogenic E. coli for 8 h.
RNA-seq on WS190, TopHat, FPKM from Partek, then DESeq.
Supplemental methods: "S. aureus-induced genes" are induced in infected vs uninfected WT at p ≤ 0.05 after multiple-test correction (825 genes).
HLH-30-dependent genes are the WT-induced genes that are not induced (p > 0.05) in hlh-30.
HLH-30-independent genes are induced at p ≤ 0.05 in both backgrounds.
637 + 188 = 825.
The lists contain induced genes only; the paper gives no list of repressed genes.

- `visvikis2014_hlh30_dependent_saureus_induced.tsv`: sheet "Table S3". Columns: FPKM means, fold change of means for WT and hlh-30, and the WT and hlh-30 p values. Unmapped: c30g4.7, c44h9.8, f45e1.3, h16o14.2, t18d3.8, y48e1c.3.
- `visvikis2014_hlh30_independent_saureus_induced.tsv`: sheet "Table S5". Same columns. Unmapped: f59a3.5, k07b1.8, t28f4.6, y111b2a.28.
- `visvikis2014_saureus_induced_wt.tsv`: sheet "Table S2". WT FPKM, fold change and p value.
- `visvikis2014_hlh30_dependent_ebox.tsv`: sheet "Table S6", "HLH-30 dependent S. aureus-induced genes containing E-box motif in 2 Kb upstream of transcription start site". Gene names only.

Caveat: the source writes names in lower case, including sequence names (for example `f21c10.10`).
The names were mapped case-insensitively against WormBase in `gene_db.sqlite`, which is newer than WS190.
Unmapped names are probably dead or renamed sequence names.
The sheet also has "Table S4" (684 genes induced in hlh-30 mutants); it is not converted.

## O'Rourke and Ruvkun 2013 (fasting)

Citation: O'Rourke EJ, Ruvkun G. MXL-3 and HLH-30 transcriptionally link lipolysis and autophagy to nutrient availability. Nat Cell Biol 2013;15(6):668-676. doi:10.1038/ncb2741.

Source: `https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fncb2741/MediaObjects/41556_2013_BFncb2741_MOESM32_ESM.xlsx` (Supplementary Table 4) saved as `raw/41556_2013_BFncb2741_MOESM32_ESM.xlsx`.

- `orourke2013_hlh30_dependent_autophagy_qpcr.tsv`: Sheet1. "HLH-30-dependent autophagic genes": genes induced after 8 h fasting of young adults in WT but not in hlh-30(tm1978), by qPCR. Columns: mean ddCT ± SEM relative to fed WT, n = 3. Six genes: atg-7, atg-9, atg-10, atg-13, atg-16.1, unc-51.

This is a qPCR candidate panel, not a genome-wide list.
The column header says hlh-30(tm1958), but the legend says tm1978.
The other supplementary tables (1 to 3, 5 to 10) are qPCR primer lists, a curated list of nutrition-responsive genes, a TF RNAi screen list, lifespans and MXL-3 ChIP-qPCR.
None of them is an HLH-30 target list.

## Lapierre et al. 2013 (autophagy and longevity)

Citation: Lapierre LR, et al. The TFEB orthologue HLH-30 regulates autophagy and modulates longevity in Caenorhabditis elegans. Nat Commun 2013;4:2267. doi:10.1038/ncomms3267. PMC3866206.

Source: `https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fncomms3267/MediaObjects/41467_2013_BFncomms3267_MOESM732_ESM.pdf` saved as `raw/41467_2013_BFncomms3267_MOESM732_ESM.pdf`.
The paper has only this PDF as supplementary material.

- `lapierre2013_clear_motif_candidates.tsv`: Supplementary Table S1, typed in by hand from the PDF. 17 C. elegans orthologues of mammalian TFEB autophagy and lysosome targets. `n_clear_sites_7of10` counts the CLEAR-like sites with at least 7 of 10 bp matching GTCACGTGAC, within 1 kb upstream to 200 bp downstream of the start codon. vps-11, vha-16 and asp-1 have 0 sites at this cutoff.

This is a candidate list chosen by orthology and checked with a motif scan.
It is not an hlh-30-dependent expression list.
The PDF prints the atg-18 cosmid as "F4E6.13"; the correct name F41E6.13 was used for mapping.

## Harvald et al. 2017 (starvation multi-omics)

Citation: Harvald EB, et al. Multi-omics analyses of starvation responses reveal a central role for lipoprotein metabolism in acute starvation survival in C. elegans. Cell Syst 2017;5(1):38-52.e4. doi:10.1016/j.cels.2017.06.004.

Source: `https://ars.els-cdn.com/content/image/1-s2.0-S2405471217302326-mmc2.xlsx` (Table S1, 40 MB) saved as `raw/1-s2.0-S2405471217302326-mmc2.xlsx`.

Design: N2 and hlh-30(tm1978), starved at mid-L4 for 0, 1, 2, 3, 4, 6 and 16 h, with label-free proteomics and RNA-seq in triplicate.
For each protein or mRNA the authors compute A, the correlation between the WT and hlh-30 time courses, and B, the integrated area between the two curves.

- `harvald2017_known_hlh30_targets.tsv`: sheets "Grove et al. HLH-30 targets_dec" (122 rows, `direction = down`) and "Grove et al. HLH-30 targets_inc" (10 rows, `direction = up`). The authors give WBGene IDs and UniProt accessions, and the B and score values for protein and mRNA where measured. 132 rows and 128 distinct genes (UniProt-to-gene duplicates, such as WBGene00004145 for pqn-62 to pqn-65).
- `harvald2017_starvation_hlh30_protein_highB_clear.tsv`: sheet "Proteins B score > 1.85 _ CLEAR", header "Proteins with B score > 1.85_with CLEAR element". These are the 1,109 high-B proteins (75th percentile) restricted to genes with a CLEAR promoter motif. Gene names only. Unmapped: R08E3.2, C18C4.5.

Uncertain: the paper text only says these were compared "with previously identified HLH-30 targets (Grove et al., 2009)".
Grove et al. 2009 ran a limma microarray of hlh-30(tm1978) vs N2 (pooled embryo, L1 and L2), with DE at log2FC ≥ 0.7 and P < 0.01.
These two sheets are most likely that DE list.
I read "dec" as lower in hlh-30, so `direction = down` means probably HLH-30-activated, but no legend confirms this.
Resolved: Grove's Table S8 does hold the microarray DE list (`new_degsets/grove2009_hlh30_microarray_de.tsv`, 123 down and 12 up). All 118 "dec" genes are in its down block and all 10 "inc" genes are in its up block, so `down` here means lower in hlh-30 (see `new_degsets/NOTES.md`).
The sheets also marked CLEAR-containing genes by cell colour, which is lost in the TSV.
The full 1,109-protein list (sheet "LOG2 ratio highscoring proteins", UniProt accessions only) is not converted, because there is no UniProt-to-WBGene mapping available locally.
B has no sign, so the CLEAR list has no direction.

## Grove et al. 2009 (bHLH network, PBM)

Citation: Grove CA, et al. A multiparameter network reveals extensive divergence between C. elegans bHLH transcription factors. Cell 2009;138(2):314-327. doi:10.1016/j.cell.2009.04.058. PMC2774807.

Source: `https://ars.els-cdn.com/content/image/1-s2.0-S0092867409005194-mmc7.xls` (Table S6, "candidate target genes") saved as `raw/1-s2.0-S0092867409005194-mmc7.xls`.
The methods PDF is `raw/1-s2.0-S0092867409005194-mmc1.pdf`.

- `grove2009_hlh30_pbm_predicted_targets.tsv`: column "HLH-30" of Table S6. 548 WBGene IDs predicted from HLH-30 protein-binding-microarray 8-mer specificity in promoters. These are predicted targets, not measured ones. 14 IDs are not in the current `gene_db.sqlite` (probably dead or merged genes), so their `gene_name` is empty.

## ENCODE / modERN HLH-30 ChIP-seq, L4

Citation: Kudron MM, et al. The ModERN resource: genome-wide binding profiles for hundreds of Drosophila and Caenorhabditis elegans transcription factors. Genetics 2018;208(3):937-949. doi:10.1534/genetics.117.300657. ENCODE experiment ENCSR246PQU.

Source: `https://www.encodeproject.org/files/ENCFF793PTZ/@@download/ENCFF793PTZ.bed.gz`, the "optimal IDR thresholded peaks" (ce11, replicates 1 and 2), saved as `raw/ENCFF793PTZ.bed.gz`.
Strain OP433 (wgIs433, HLH-30::eGFP fosmid transgene), L4 hermaphrodites, anti-GFP ChIP.
The file has 7,168 peaks.

- `encode_modern_hlh30_l4_chip_promoter_bound.tsv`: built by me, not by the authors. A protein-coding gene counts as bound when a peak summit lies from 1,000 bp upstream to 200 bp downstream of any of its transcript starts in `work/*/Caenorhabditis_elegans.WBcel235.115.gtf` (strand-aware). `summit_to_tss` is the nearest such distance (negative means upstream). `signal` is the narrowPeak signalValue.

Caveat: one peak can be assigned to two genes (bidirectional promoters, operons), so this list over-counts real targets.
The window is my own choice, not one taken from a paper.
Condition: fed L4 with a multicopy transgene, not stress-activated HLH-30.

## Sources tried and not obtained

- Lin XX, et al. 2018, Nat Commun 9:4400, doi:10.1038/s41467-018-06624-0, PMC6199276. The only supplementary files are `41467_2018_6624_MOESM1_ESM.pdf` (16 pages: strain list and lifespan statistics, Supplementary Tables 1 and 2) and the peer review file. The paper publishes no gene list for HLH-30 ChIP peaks or for hlh-30-dependent DEGs. The data availability statement gives only raw reads (SRA SRP152334 and SRP017927) and data "upon reasonable request". To get Lin's lists, the reads would need to be reprocessed.
- Visvikis 2014 through PMC and Europe PMC: the PMC `/bin/` links return a JavaScript proof-of-work page to curl, Europe PMC returns HTTP 403, and the Europe PMC supplementaryFiles API says "not open access". The files came from the Elsevier CDN instead (see above).
- Sohn J, et al. 2025, Autophagy, "HLH-30/TFEB mediates sexual dimorphism in immunity" (doi:10.1080/15548627.2024.2375779; earlier listed here as "Goswamy 2025"): not open access in Europe PMC, not in the PMC open-data bucket. Later obtained from the publisher's figshare; see `new_degsets/NOTES.md`.
- Dall KB, et al. 2021, Aging Cell 20:e13342 (PMC8045935): supplements are metabolomics, lipidomics and survival tables only. Its transcriptome is Harvald 2017, used above.
