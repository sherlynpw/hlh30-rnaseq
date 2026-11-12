# Autophagy and lysosome genes linked to HLH-30

Three gene lists:

1. Autophagy and lysosome genes that published papers link to HLH-30. These go beyond the two lists in `../` (O'Rourke 2013 qPCR, Lapierre 2013 CLEAR candidates).
2. A curated list of core C. elegans autophagy genes, taken from a review table. It does not depend on HLH-30.
3. The GO:0006914 (autophagy) gene set from org.Ce.eg.db.

Downloads are in `raw/`.
`build_lists.py` turns them into the tidy TSVs in this folder:

```
Rscript export_go.R
uv run --with pandas python build_lists.py
```

`export_go.R` writes `raw/go0006914_orgCeegdb.tsv`, and `build_lists.py` reads that file.
Names map to WBGene IDs the same way as `../build_lists.py`: a case-insensitive match against `gene_name`, then `sequence_name`, then `other_names` in `gene_db.sqlite`.
Some author spellings need a fixed target. They are listed in `ALIASES` in `build_lists.py`, and those rows get `map_method = alias`.

## Summary

| file | meaning | rows | mapped (distinct IDs) |
|---|---|---|---|
| hlh30_autophagy_evidence.tsv | one row per gene, paper, condition and evidence kind | 208 | 207 (65) |
| hlh30_autophagy_genes_summary.tsv | one row per gene: papers per evidence kind | 65 | 64 (64) |
| core_autophagy_zhang_melendez2025.tsv | Zhang & Meléndez 2025 Table 1 | 71 | 71 (69) |
| go0006914_autophagy_orgCeegdb.tsv | GO:0006914 with descendants (GOALL), org.Ce.eg.db 3.22.0 | 127 | 127 (127) |

## Evidence list

`evidence_curated.tsv` is the input.
I built it by hand from the papers.
Each row cites a figure or table panel and has a short verbatim `quote` from the legend or the results.
Where gene names appear only on a figure axis, I read them from the figure image (publisher JPG, or a rendered supplementary PDF page), and the `notes` column says so.
`build_lists.py` adds `gene_id`, `gene_name`, `gene_class`, `in_core_list`, `in_go0006914` and `map_method` and writes `hlh30_autophagy_evidence.tsv`.

Columns:

- `evidence` holds one of five kinds:
  - `expr_reduced_hlh30`: the mRNA or a transcriptional reporter is lower, or its induction is lost or blunted, in an hlh-30 mutant or after hlh-30 RNAi. The `result` column says which.
  - `induced_hlh30_oe`: higher with HLH-30 overexpression.
  - `chip_binding`: HLH-30 ChIP at the promoter, shown in that paper.
  - `motif_only`: a predicted CLEAR, E-box or HLH-30 site, with no binding data. This is the weakest kind.
  - `not_dependent`: the paper tested hlh-30 loss and the change was still there. Kept as a negative.
- `hlh30_perturbation`, `condition`, `stage_tissue`, `assay`, `result`, `figure`, `paper`, `pmcid`, `doi`, `quote` and `notes` describe each claim.
- `gene_class` is the `role_class` from the core list (`core_machinery`, `maturation_fusion`, `selective_autophagy`, ...), or `lysosome_or_other` when the gene is not in the core list.

The summary file counts papers per evidence kind; `not_dependent` papers are listed separately.
`strongest_evidence` orders the kinds as chip_binding, then expr_reduced_hlh30, then induced_hlh30_oe, then motif_only.
`n_papers_functional` counts papers with any kind except `motif_only`.

What counts: only claims a paper shows in its own figure or table, with hlh-30 perturbed (mutant, RNAi or overexpression).

What does not count:
- an expression change in another mutant or treatment where hlh-30 was not perturbed
- GFP::LGG-1 puncta, tandem reporters or translational fusions (these measure flux, not transcription)
- statements that cite earlier work

Distinct genes per evidence kind (mapped IDs):

| evidence | genes | of which in the core autophagy list |
|---|---|---|
| expr_reduced_hlh30 | 56 | 21 |
| induced_hlh30_oe | 14 | 4 |
| chip_binding | 2 | 0 |
| motif_only | 22 | 8 |
| not_dependent (at least one paper) | 12 | |

### Papers used

All papers used hlh-30(tm1978) unless the table says otherwise.

| paper | condition | panels | rows |
|---|---|---|---|
| Lapierre 2013 Nat Commun, doi:10.1038/ncomms3267 | glp-1 (Fig. 1c), tor RNAi (Fig. 2b,c), glp-1 + adult hlh-30 RNAi (Supp. Fig. S1a), WT + hlh-30 RNAi (S1b), basal tm1978 (S1c), HLH-30::GFP overexpression (S6b), CLEAR sites (Table S1) | Fig. 1c, 2b-c, S1a-c, S6b, Table S1 | 53 |
| O'Rourke 2013 Nat Cell Biol, doi:10.1038/ncb2741 | fasting (Fig. 5b, lgg-1/lgg-2/atg-16.2; Supp. Table 4, 8 h), lipl genes after 5 h fasting (Fig. 4a), mxl-3 background (Fig. 4g), HLH-30::GFP ChIP-qPCR at the lipl-1 and lipl-3 promoters after 6 h fasting (Fig. 4f), overexpression (Supp. Fig. S6c) | | 23 |
| Visvikis 2014 Immunity, doi:10.1016/j.immuni.2014.05.002 | S. aureus 8 h; rescue with HLH-30::GFP | Fig. 5A-C, S4A, S4H | 17 |
| Nakamura 2016 Nat Commun, doi:10.1038/ncomms10944 | glp-1;hlh-30 vs glp-1 (Fig. 6a,c), let-363 RNAi (Supp. Fig. 7a) | | 6 |
| Kumsta 2017 Nat Commun, doi:10.1038/ncomms14337 | hormetic heat shock, 1 h at 36 C | Supp. Fig. 9a | 4 |
| Chen HD 2017 Autophagy, doi:10.1080/15548627.2016.1256933 | Cry5B pore-forming toxin: qRT-PCR 3 h (Fig. 6G, 7D); RNA-seq 1 h (Fig. 7A-B, Table S1); predicted sites (Fig. 6A) | | 16 |
| Bennett 2017 PLoS Genet, doi:10.1371/journal.pgen.1006695 | tald-1 RNAi | Fig. 6E-F | 4 |
| Franco-Juárez 2018 Aging, doi:10.18632/aging.101577 | high-glucose diet | Fig. 2D | 4 |
| Murphy 2019 PLoS Biol, doi:10.1371/journal.pbio.3000245 | L1 starvation 33 h, then refeeding on OP50 or CeMM | Fig. 4A, 5C-N, S4, S18, S20 | 29 |
| Imanikia 2019 Curr Biol, doi:10.1016/j.cub.2019.06.031 | neuronal xbp-1s, intestine-only hlh-30 RNAi; CLEAR sites | Fig. S5C, S3A | 13 |
| Dinić 2021 Aging, doi:10.18632/aging.202885 | heat-killed Lb. fermentum BGHV110, 6 h | Fig. 3C, 4C | 8 |
| Kumar 2022 Sci Adv, doi:10.1126/sciadv.abj1604 | xpo-1 RNAi | Fig. S3F | 1 |
| Possik 2023 Nat Commun, doi:10.1038/s41467-023-40857-y | pgph-2 overexpression, RNA-seq | Supp. Fig. 4a, Supp. Data 10 | 3 |
| Chamoli 2023 Nat Aging, doi:10.1038/s43587-023-00524-9 | mitophagy-inducing coumarin (MIC) | Ext. Data Fig. 2d | 6 |
| Metcalf 2024 iScience, doi:10.1016/j.isci.2024.109354 | glial XBP-1s, RNA-seq | Fig. 3G, S3G, Table S5 | 8 |
| Lin H 2024 Front Pharmacol, doi:10.3389/fphar.2024.1433030 | hlh-30 RNAi (Fig. 3A); HLH-30 overexpression in the Abeta strain GMC101 (Supp. Fig. S1) | | 7 |
| Zhong 2025 PLoS Biol, doi:10.1371/journal.pbio.3002957 | rab-7 transcriptional reporter in PVD neurons; predicted sites | Fig. 3C-F | 2 |
| Li TY 2025 Nat Cell Biol, doi:10.1038/s41556-025-01693-y | vha-6 RNAi lysosomal surveillance response; cpr-5 and cpr-8 are not hlh-30-dependent | Fig. 3c, Ext. Data Fig. 3d, 10a | 4 (all not_dependent) |

Caveats, also in the `notes` column:

- Some significance calls were read from bar-chart stars in figure images, because the text does not state them. This applies to the Lapierre not_dependent rows, Murphy S4/S20, Kumsta S9a, Chamoli ED2d and Dinić 3C.
- Lapierre `ctsa*` is cosmid C08H9.1 (Table S1), now named ctsa-3.2. The old name ctsa-1 is an alias of ctsa-4.1, so `ctsa*` is mapped by cosmid.
- Chamoli writes only "ctsa" and gives no cosmid, so that row is unmapped.
- Metcalf's ctsa-2 was mapped by the current name (K10C2.1). Y40D12A.2 also has ctsa-2 as an old alias.
- Chen 2017 Table S1: the "537 hlh-30-dependent Cry5B response genes" are defined in the text as up with Cry5B in N2 and down in hlh-30. The log2FC column is not labeled with a direction. I took the autophagy genes in it (lgg-1, lgg-2, atg-2, atg-3, atg-9, atg-18, epg-6) as `expr_reduced_hlh30`.
- O'Rourke Supp. Table 4 prints the allele as tm1958, but the legend says tm1978.
- Bennett sqst-1 is still about 3.5-fold induced in hlh-30 tald-1(RNAi), with no significance mark. The authors still call it HLH-30-dependent.
- Dinić lipl-1 is already about 2-fold higher in untreated hlh-30. The lost induction may be a ceiling effect.
- Dinić atg-18 goes down in hlh-30 after BGHV110, rather than staying flat.
- Nakamura: unc-51 is higher in glp-1;hlh-30 than in glp-1 (Fig. 6a; Supp. Fig. 6e). No evidence kind fits this, so it is not in the TSV.
- Nakamura Supplementary Data 1-2 are DAVID enrichment tables with an unclear sign, so no rows came from them.
- Murphy vha-8, vha-9 and vha-14 (starvation rows) compare starved with fed hlh-30, not WT with hlh-30.
- Possik's 5 lysozyme genes (ilys-2, lys-4 to lys-7) were dropped as out of scope.
- Metcalf: cdr-4, asp-6, asp-17 and ctsa-2 are lysosomal by the paper's annotation, but they are not on its 36-gene HLH-30 target list.
- Imanikia: cdr-1 is lysosomal by the authors' grouping.

### Papers checked with no qualifying claim

- Silvestrini 2018 Cell Rep: the qPCR of lgg-1, lgg-2 and sqst-1 after xpo-1 RNAi was done in WT only. hlh-30 was tested only with puncta and lifespan.
- Lin XX 2018 Nat Commun: no autophagy or lysosome gene is named in the ChIP or expression figures, only GO terms.
- Seah 2016 Autophagy: no expression was measured in hlh-30.
- Ji 2020 PLoS Genet: atg-9, sqst-1 and the lipl genes were measured in alfa-1 only.
- Zhou 2019 Cell Rep: no expression was measured in hlh-30.
- Chen J 2017 eLife (metformin): no expression was measured in hlh-30.
- Wang C 2017 PNAS: the TFEB target genes were measured in mammalian cells only.
- Shioda 2023: the gene expression was done with mml-1/mxl-2 RNAi only.
- Schiffer 2021: the vha genes were measured in mir-1 mutants only.
- Ghaddar 2023: the only hlh-30-dependent gene shown is adh-1.
- Nonninger 2025: no expression of autophagy or lysosome genes was measured in hlh-30.
- Shalash 2025: the hlh-30-dependent genes shown are chaperones only.
- The 2026 L1 quiescence chromatin paper (PMC13646826): no autophagy or lysosome genes are named.
- Kumar 2023 (lipid droplets): hlh-30 was used for lifespan only.
- Taylor 2024 (MANF): hlh-30 was not perturbed.
- Ricaurte-Perez 2025: "vps-34 is under transcriptional control by HLH-30" is not backed by a panel with hlh-30 perturbed.
- Tan 2025 (AA-20): WT only.
- Hsieh 2024 (PRMT-7): only prmt-7 was perturbed.
- Lee 2025 (pelo-1): only GSEA, with no named genes.
- Lazaro-Pena 2023: puncta only.

## Search

- Europe PMC query, saved in `raw/europepmc_search_hlh30_autophagy.json` (499 hits, 459 with a PMCID): `("hlh-30" OR "HLH-30") AND (autophagy OR autophagic OR "lgg-1" OR "lgg-2" OR "bec-1" OR "sqst-1" OR "atg-18" OR "unc-51" OR lysosome OR lysosomal)`.
- Full text: Europe PMC `fullTextXML` (335 papers). For the other 124 I used NCBI efetch (`db=pmc`), which returned the full body for author manuscripts such as Lapierre 2013, O'Rourke 2013, Visvikis 2014 and Silvestrini 2018. About 58 returned only front matter.
- I screened every sentence and figure legend that names hlh-30 or TFEB, an autophagy or lysosome gene and an expression or ChIP word. I then read the candidate papers.
- I also screened the papers citing Lapierre 2013 or O'Rourke 2013 that were not in the search (`raw/europepmc_citing_lapierre2013_orourke2013.json`; 160 with full text). They added no new claims.
- `raw/europepmc_xml/` and `raw/pmc_efetch/` keep only the XML of the papers listed above. `raw/<PMCID>/` holds the supplementary files and figure images each row was read from.

## Core autophagy list

Citation: Zhang H, Meléndez A. Conserved components of the macroautophagy machinery in Caenorhabditis elegans. Genetics 2025;229(4):iyaf007. doi:10.1093/genetics/iyaf007. PMC12005284.

- `core_autophagy_zhang_melendez2025.tsv` is Table 1, "Autophagy-related genes and autophagy regulators in yeast, C. elegans, and mammals", parsed from `raw/europepmc_xml/PMC12005284.xml`.
- `category` is the table's own group, such as "Atg1/ULK complex" or "HOPS complex, for vesicle fusion".
- `role_class` is my grouping of those categories: `core_machinery`, `maturation_fusion`, `selective_autophagy`, `autophagy_related_or_regulator` or `transcriptional_regulator`.
- atg-7 appears twice because it is in both conjugation systems.
- 71 rows give 69 genes, and 23 genes are `core_machinery`.
- `in_zhang2015_table1` marks the genes that are also in Table 1 of Zhang H, et al. Guidelines for monitoring autophagy in C. elegans. Autophagy 2015;11(1):9-27. doi:10.1080/15548627.2014.1003478. That table has the 26 core ATG and EPG genes and no regulators.
- Names in the table that needed a fix:
  - `epg-1` is now atg-13 (other_names).
  - `lpla-2` is now plag-15 (other_names).
  - `daf-21` is now hsp-90 (other_names).
  - `Y34BA.2` is a typo for Y34B4A.2.
  - `impk-1` is ipmk-1.
  - `fundc-1` is fndc-1.
  - `rub-1` was mapped to Y56A3A.16, which WormBase calls the RUBCN ortholog.
  - `allo-1a` and `allo-1b` are both allo-1.

## GO:0006914

- `go0006914_autophagy_orgCeegdb.tsv` comes from `AnnotationDbi::select(org.Ce.eg.db, keys = "GO:0006914", keytype = "GOALL")`. It includes descendant terms such as mitophagy and the autophagosome assembly terms.
- Source: org.Ce.eg.db 3.22.0, GO source date 2025-07-22.
- It has 127 genes. 32 of them have IEA evidence only (`iea_only`).
- 47 of the 69 core genes are in the GO set.
- These core genes are missing from the GO set: rab-7, syx-17, vps-11, vps-16, vps-33.1, vps-33.2, snap-29, vamp-7, ogt-1, plag-15, hsp-70, hsp-90, rheb-1, rict-1, daf-16, rab-11.1, rab-11.2, vps-4, ipmk-1, T23G11.7, Y34B4A.2 and Y56A3A.16.

## Not obtained

- Gerisch B, et al. 2020 Dev Cell (HLH-30 and reproductive quiescence): not in PMC.
- Sohn J, et al. 2025 Autophagy (sexual dimorphism, PMC11759534): not open access in Europe PMC; not read.
- These 2026 papers have no PMCID yet:
  - Hsieh 2026 (methyl-HLH-30)
  - Cubillas 2026 (zinc and lysosomes)
  - Schmauck-Medina 2026 (ATG-18 and HLH-30)
  - Meng 2026 (puerarin)
- Franco-Juárez 2018: aging-us.com blocks curl with a Cloudflare challenge. The Fig. 2D gene set comes from the text, and the per-gene bars were not checked against the image.
- Silvestrini 2018 supplement: the Europe PMC zip had one empty image, and the NIHMS supplement link returned HTML.
- Gelino 2016 PLoS Genet did not match the search. I did not check it further.
