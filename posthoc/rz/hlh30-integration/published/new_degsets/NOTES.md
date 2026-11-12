# Additional published HLH-30 gene lists

This folder holds lists found after the first set in `../README.md`.
Original files are in `raw/`.
`build_new.py` turns them into the tidy TSVs in this folder:

```
uv run --with pandas --with openpyxl --with xlrd --with python-docx python build_new.py
```

## TSV format

The columns are the same as in `../README.md`: `gene_id`, `gene_name`, `source_gene`, `direction`, extra score columns, `map_method`.
When the authors give a WBGene ID, `map_method = source` and `gene_name` comes from `gene_db.sqlite`.
If that ID is no longer in `gene_db.sqlite` (dead or merged gene), `gene_name` is empty.
Otherwise names are mapped as in `build_lists.py` (`gene_name`, `sequence_name`, `other_names`, case-insensitive).

`direction` has one of two meanings, set per file below:

- Mutant vs control lists (Grove, Sohn, Gerisch, Shioda, Lin, Munoz-Barrera, Metcalf): `up` and `down` mean higher or lower in the hlh-30 loss-of-function sample than in its control.
- Stimulus-induced, hlh-30-dependent lists (Chen, Wani, Hsieh, Possik, Murphy group D): `up` means induced in the control background by the stimulus, and not induced (or less induced) without HLH-30. This is the meaning `up` has in the Visvikis lists of `../README.md`.

## Summary

Rows, mapped rows and distinct IDs are from the last run of `build_new.py`.

| file | meaning | rows | mapped (distinct) | direction counts |
|---|---|---|---|---|
| grove2009_hlh30_microarray_de.tsv | hlh-30 vs N2 microarray, mixed embryo/L1/L2 | 135 | 135 (135) | down 123, up 12 |
| sohn2025_hlh30_down_in_males.tsv | male-enriched genes lower in hlh-30 males | 322 | 322 (322) | down 322 |
| chen2017_hlh30_dependent_cry5b_induced.tsv | Cry5B-induced in WT, hlh-30-dependent | 521 | 521 (521) | up 521 |
| wani2021_hlh30_dependent_saureus_induced.tsv | S. aureus-induced (vs starved) in WT only | 99 | 99 (99) | up 99 |
| wani2021_hlh30_independent_saureus_induced.tsv | S. aureus-induced in WT and hlh-30 | 6 | 6 (6) | up 6 |
| hsieh2026_hlh30_r38_dependent_cry5b_induced.tsv | Cry5B-induced, lost with HLH-30(R38K) rescue | 314 | 314 (314) | up 314 |
| gerisch2020_hlh30_vs_wt_ard_de.tsv | hlh-30 vs N2, adult reproductive diapause (starved adults) | 3666 | 3666 (3666) | down 1884, up 1782 |
| gerisch2020_hlh30_vs_wt_al_de.tsv | hlh-30 vs N2, fed adults | 124 | 124 (124) | down 111, up 13 |
| possik2023_hlh30_dependent_pgph2_induced.tsv | induced by pgph-2 overexpression, hlh-30-dependent | 214 | 213 (213) | up 214 |
| shioda2023_glp1_neuronal_hlh30_rnai_de.tsv | neuronal hlh-30 RNAi vs control in glp-1 | 93 | 93 (93) | down 70, up 23 |
| lin2024_gmc101_hlh30_ko_de.tsv | hlh-30 KO vs control in Abeta strain GMC101 | 62 | 60 (60) | up 49, down 13 |
| nonninger2025_hlh30_vs_wt_ard_log2fc_all_genes.tsv | all genes with ARD log2FC hlh-30/N2 (no p values) | 17694 | 17579 (17572) | NA |
| munozbarrera2026_hlh30_vs_wt_starved_l1_de.tsv | hlh-30 vs N2, starved L1 | 7899 | 7899 (7899) | down 3996, up 3903 |
| murphy2019_hlh30_starved_l1_kegg_clusters.tsv | heatmap clusters, starved L1, KEGG-annotated members only | 185 | 185 (185) | down 78, up 53, up_in_hlh30 54 |
| metcalf2024_glial_xbp1s_hlh30_vs_xbp1s_de.tsv | glial XBP-1s;hlh-30 vs glial XBP-1s, day 2 adults | 37 | 37 (37) | up 24, down 13 |
| encode_modern_hlh30_late_embryo_chip_promoter_bound.tsv | HLH-30::GFP ChIP summit at promoter, late embryo | 1328 | 1328 (1328) | bound 1328 |

## Grove et al. 2009 (microarray, Table S8)

Citation: Grove CA, et al. A multiparameter network reveals extensive divergence between C. elegans bHLH transcription factors. Cell 2009;138(2):314-327. doi:10.1016/j.cell.2009.04.058. GEO GSE15762.

Source: `https://ars.els-cdn.com/content/image/1-s2.0-S0092867409005194-mmc9.xls` (Table S8, "HLH-30 target genes"), saved as `raw/1-s2.0-S0092867409005194-mmc9.xls` (36 KB).

Design: N2 vs hlh-30(tm1978), mixed-stage populations (see `../README.md`: pooled embryo, L1 and L2), two-colour microarray, limma.
The README gives the paper's DE call as log2FC ≥ 0.7 and P < 0.01; the table holds P.Value and adj.P.Val, and the smallest |logFC| in it is 0.667.
Almost no gene passes an adjusted P threshold, so this is a nominal-P list.

- `grove2009_hlh30_microarray_de.tsv`: block "Genes with increased expression in hlh-30(tm1978) mutant" gives `direction = up` (12), block "Genes with decreased expression" gives `down` (123). Columns `logFC`, `pvalue`, `padj` as given. WBGene00018464 (F45E1.3) is not in `gene_db.sqlite`.

This resolves the uncertainty in `../README.md` about the Harvald "Grove et al. HLH-30 targets" sheets: all 118 "dec" genes are in the `down` block here and all 10 "inc" genes are in the `up` block.
So Harvald's "dec" means lower in hlh-30.
`../README.md` says Grove's supplement has no microarray DE list; Table S8 (mmc9) is that list.

## Sohn et al. 2025 (males, sexual dimorphism in immunity)

Citation: Sohn J, et al. HLH-30/TFEB mediates sexual dimorphism in immunity in Caenorhabditis elegans. Autophagy 2025. doi:10.1080/15548627.2024.2375779. PMC11759534. GEO GSE245430.
`../README.md` lists this paper as "Goswamy 2025"; the first author is Sohn (Lee SJV lab).

Source: Taylor & Francis figshare (doi:10.6084/m9.figshare.26176576), `https://ndownloader.figshare.com/files/47426455`, Table S3 "Dataset - RNA seq analysis", saved as `raw/kaup_a_2375779_sm6158.xlsx` (108 KB).

Design: N2 and hlh-30(tm1978), hermaphrodites and males, fed on E. coli, n = 3 per group. RSEM, upper-quartile + RUVg normalisation, DESeq2 Wald with BH.

- `sohn2025_hlh30_down_in_males.tsv`: sheet "Fig. 5A", "322 genes downregulated by hlh-30(-) in males (absolute fold change > 1.5, padj < 0.05)". Every row is also male-enriched (WT male vs WT hermaphrodite log2FC > 0.585, padj < 0.05). `direction = down` (lower in hlh-30 males). Columns: hlh-30 male vs WT male log2FC and padj, WT male vs hermaphrodite log2FC and padj.

Caveat: the list is restricted to male-enriched genes. There is no hermaphrodite hlh-30 vs WT list; GEO has only per-sample RSEM counts.

## Chen et al. 2017 (Cry5B pore-forming toxin)

Citation: Chen HD, et al. HLH-30/TFEB-mediated autophagy functions in a cell-autonomous manner for epithelium intrinsic cellular defense against bacterial pore-forming toxin in C. elegans. Autophagy 2017;13(2):371-385. doi:10.1080/15548627.2016.1256933. PMC5324838. GEO GSE78878.

Source: figshare doi:10.6084/m9.figshare.4248818, `https://ndownloader.figshare.com/files/6927929`, Table S1 "The list of hlh-30-dependent Cry5B response genes", saved as `raw/2016AUTO0201R2-s03.xlsx` (51 KB).
GEO processed files for the four pairwise comparisons are also saved (`raw/GSE78878_*_score.txt.gz`, 0.36-0.39 MB each; keyed by RefSeq NM_ IDs; not converted).

Design: N2 and hlh-30(tm1978) L4 to young adults on E. coli expressing Cry5B or the empty vector (pQE9). Exposure 1 h in the text, 3 h in GEO. Illumina MiSeq, apparently one library per condition (GEO has 4 samples).
Criterion: up in N2 Cry5B vs N2 control, and down in hlh-30 Cry5B vs N2 Cry5B, and down in hlh-30 control vs N2 control, each at P < 0.05 without FDR correction.

- `chen2017_hlh30_dependent_cry5b_induced.tsv`: 544 source rows, 521 distinct WBGene IDs (duplicate transcripts; the row with the smallest p is kept). The paper text says 537 genes. `direction = up`. `log2fc_normalized` is the column as given; it is negative for all but 2 rows, so it is probably hlh-30 vs N2 (the legend does not say). 3 IDs are not in `gene_db.sqlite`.

Caveat: n = 1 per condition and nominal P values.

## Wani et al. 2021 (S. aureus, NHR-49 and HLH-30)

Citation: Wani KA, Goswamy D, et al. NHR-49/PPAR-α and HLH-30/TFEB cooperate for C. elegans host defense via a flavin-containing monooxygenase. eLife 2021;10:e62775. doi:10.7554/eLife.62775. PMC8139828. SRA PRJNA727590, PRJNA728096.

Source: `https://cdn.elifesciences.org/articles/62775/elife-62775-supp1-v3.xlsx` (0.93 MB) and `https://cdn.elifesciences.org/articles/62775/elife-62775-supp3-v3.xlsx` (0.49 MB), saved under the same names in `raw/`.

Design: WT and hlh-30(tm1978) young adults, S. aureus SH1000 vs starvation (the reference condition), n = 4, DESeq2. WT DEGs appear to use padj < 0.01 (largest padj in the sheet is 0.00999).

- `wani2021_hlh30_dependent_saureus_induced.tsv`: supp3, sheet "Overlap_UP_WT_v_hlh-30(-)", column "Wild type ONLY" (99 genes): up in WT on S. aureus vs starved, not up in hlh-30. `wt_log2fc` and `wt_padj` are looked up by name in supp1 sheet "up+down_WT_Starv_v_Saur".
- `wani2021_hlh30_independent_saureus_induced.tsv`: column "Wild type & hlh-30(-)" (6 genes).

Caveat: the reference is starved worms, not worms on E. coli, so this is infection vs starvation.

## Hsieh et al. 2026 (methyl-HLH-30, Cry5B)

Citation: Hsieh HC, et al. Methyl-HLH-30/TFEB orchestrates holistic intestinal recuperation against bacterial pore-forming toxin intoxication. Autophagy 2026. doi:10.1080/15548627.2026.2735736. ArrayExpress E-MTAB-15698.

Source: figshare doi:10.6084/m9.figshare.33872681, `https://ndownloader.figshare.com/files/68765533`, saved as `raw/kaup_a_2735736_sm2037.docx` (24.6 MB). Table S1, "The methyl-HLH-30 R38-dependent PFT-induced genes", is the first table in the docx.

Design (from E-MTAB-15698): hlh-30(tm1978) rescued with intestinal HLH-30 (app-1p::hlh-30), WT vs R38K (methylation-site mutant), L4 to young adult, Cry5B vs pQE9 for 3 h, n = 3, NovaSeq, Salmon + DESeq2. The thresholds are in the paywalled main text and were not read.

- `hsieh2026_hlh30_r38_dependent_cry5b_induced.tsv`: 314 genes with the authors' row number in `rank`. `direction = up`.

Caveat: this is a point mutant vs WT rescue, not an hlh-30 null vs WT.

## Gerisch et al. 2020 (adult reproductive diapause)

Citation: Gerisch B, et al. HLH-30/TFEB is a master regulator of reproductive quiescence. Dev Cell 2020;53(3):316-329. doi:10.1016/j.devcel.2020.03.014. GEO GSE126632.

Source: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126632/suppl/GSE126632_gene.expression.tsv.gz`, saved as `raw/GSE126632_gene.expression.tsv.gz` (5.4 MB). The Elsevier supplement holds only lifespans and GO tables.

Design: N2 vs hlh-30(tm1978), fed adults (AL) and 24 h adult reproductive diapause (ARD, starved from mid-L4/young adult), RNA-seq, n = 4. The file is the authors' Cuffdiff gene_exp.diff with all six pairwise comparisons.
`log2fc` is log2(hlh-30 / N2) (checked: `value_2 / value_1` with sample_2 = h30).

- `gerisch2020_hlh30_vs_wt_ard_de.tsv`: comparison N2_ARD vs h30_ARD, rows with Cuffdiff `significant = yes` (q < 0.05).
- `gerisch2020_hlh30_vs_wt_al_de.tsv`: comparison N2_AL vs h30_AL, same rule.

`paper_go_cutoff` is TRUE when |log2FC| ≥ 2 and p < 0.05, the cutoff the Table S3 legend gives for the authors' GO analysis. Within the q < 0.05 rows this is 1,078 genes for ARD (564 down, 514 up) and 35 for AL (28 down, 7 up).
The list uses Cuffdiff's own call; the paper's text does not state a DEG list.
44 ARD and 2 AL IDs are not in `gene_db.sqlite`.

## Possik et al. 2023 (PGPH-2 overexpression)

Citation: Possik E, et al. Glycerol 3-phosphate phosphatase/PGPH-2 counters metabolic stress and promotes healthy aging via a glycogen sensing-AMPK-HLH-30-autophagy axis in C. elegans. Nat Commun 2023;14:5214. doi:10.1038/s41467-023-40857-y. PMC10457390. GEO GSE228784.

Source: Supplementary Table S8, `41467_2023_40857_MOESM11_ESM.xlsx` (25 KB), from static-content.springer.com. Table S2 (`41467_2023_40857_MOESM5_ESM.xlsx`, 39 KB, the 349-gene parent set) is saved but not converted. GEO raw counts `GSE228784_gene_count.txt.gz` (1.7 MB) are saved but not converted.

Design: N2, pgph-2 overexpression, pgph-2 o/e;hlh-30(tm1978), hlh-30(tm1978); L4/young adults on NGM or 2% glucose; DESeq2, n = 4.
The 214 genes are those of the 349 genes up in pgph-2 o/e vs WT on NGM (FC > 2, P < 0.05 in the legend; padj ≤ 0.05 in the methods) that go down in pgph-2 o/e;hlh-30 vs pgph-2 o/e.

- `possik2023_hlh30_dependent_pgph2_induced.tsv`: no scores. `direction = up`. Y71F9B.14 is unmapped.

Caveat: there is no author DE table for plain hlh-30 vs WT, although the counts include those samples.

## Shioda et al. 2023 (neuronal hlh-30 RNAi, glp-1)

Citation: Shioda T, et al. Neuronal MML-1/MXL-2 regulates systemic aging via glutamate transporter and cell nonautonomous autophagic and peroxidase activity. PNAS 2023;120(39):e2221553120. doi:10.1073/pnas.2221553120. PMC10523562. GEO GSE234429.

Source: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234429/suppl/GSE234429_gTU_luc_vs_gTU_hlh-30i.xlsx` (3.7 MB). Same as sheet "RPKM_neuronal hlh-30" of PNAS Dataset S01. `GSE234429_TU_luc_vs_gTU_luc.xlsx` (3.8 MB, glp-1 vs control, no hlh-30 arm) is also saved, not converted.

Design: glp-1(e2141) with neuron-specific RNAi (TU3401), hlh-30 RNAi vs luciferase RNAi from adult day 1, n = 3. The table gives RPKM per sample, logFC, p.value and FDR; the DE method was not checked.

- `shioda2023_glp1_neuronal_hlh30_rnai_de.tsv`: FDR < 0.05. `up`/`down` = higher/lower after neuronal hlh-30 RNAi.

Caveat: the FDR < 0.05 cutoff is mine; the paper text was not checked for the authors' cutoff. Knockdown is neuron-only, in a germline-less background.

## Lin H et al. 2024 (Abeta strain GMC101)

Citation: Lin H, et al. HLH-30/TFEB modulates autophagy to improve proteostasis in Aβ transgenic Caenorhabditis elegans. Front Pharmacol 2024;15:1433030. doi:10.3389/fphar.2024.1433030. PMC11392864. GEO GSE272146.

Source: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/GSE272146/suppl/GSE272146_GMC101mutvsGMC101_deg.txt.gz` (1.8 MB). `GSE272146_GMC101oxvsGMC101_deg.txt.gz` (1.8 MB, hlh-30 overexpression) is saved, not converted.

Design: hlh-30(syb6883) CRISPR knockout vs control, both in GMC101 (muscle Abeta1-42) after Abeta induction, n = 2, DESeq2.
The paper's threshold is "twofold change and FDR < 0.05".

- `lin2024_gmc101_hlh30_ko_de.tsv`: padj < 0.05 and |log2FC| ≥ 1. `up`/`down` = higher/lower in the knockout. Unmapped: novel.135, novel.69 (novel loci).

Caveat: n = 2 and an Abeta transgenic background.

## Nonninger et al. 2025 (diapause, TFEB-TGF-beta axis)

Citation: Nonninger TJ, et al. A TFEB-TGFβ axis systemically regulates diapause, stem cell resilience and protects against a senescence-like state. Nat Aging 2025. doi:10.1038/s43587-025-00911-4. PMC12270908. GEO GSE291659.

Source: Extended Data Fig. 5 source data, `43587_2025_911_MOESM12_ESM.xlsx` (0.69 MB). GEO `GSE291659_masterTable_annotated.tsv.gz` (4.9 MB, normalized counts only) is saved, not converted.

Design: N2, hlh-30(tm1978), daf-1(m40), daf-1;hlh-30 at 48 h adult reproductive diapause and after refeeding, n = 3-4, DESeq2 (padj < 0.05 in the paper).

- `nonninger2025_hlh30_vs_wt_ard_log2fc_all_genes.tsv`: sheet "Extended Figure 5a", all 17,694 genes with `ard_log2fc_hlh30_vs_n2` and `ard_log2fc_hlh30daf1_vs_hlh30`. `direction = NA`. Use it as a ranked list.

Caveat: no per-gene p values are published, so this cannot be thresholded the authors' way. 115 names are unmapped.

## Munoz-Barrera et al. 2026 (starved L1, chromatin)

Citation: Muñoz-Barrera M, et al. HLH-30/TFEB is necessary for chromatin reorganization and maintenance of cell quiescence during starvation in C. elegans. Nucleic Acids Res 2026. PMC13646826. bioRxiv doi:10.1101/2025.10.31.685810. GEO GSE313873.

Source: PMC open-access bucket `https://pmc-oa-opendata.s3.amazonaws.com/PMC13646826.1/gkag949_supplemental_files.zip` (18 MB). File S2 from that zip is saved as `raw/gkag949_File_S2_RNAseq_Rev.xlsx` (10.7 MB).

Design: N2, daf-16(mu86), hlh-30(tm1978) and daf-16;hlh-30 L1 larvae on day 1 of starvation in M9, n = 3, NextSeq 500, STAR, featureCounts, DESeq2.
log2FC is hlh-30 / WT (hlh-30 itself has log2FC -0.84).

- `munozbarrera2026_hlh30_vs_wt_starved_l1_de.tsv`: sheet "deseq2_hlh-30_wt", padj < 0.05 (the authors' DEG call, no fold-change cutoff). 4,032 of these (2,069 up, 1,963 down) have |log2FC| ≥ 1 if a stricter list is wanted.

The sheets for daf-16 and daf-16;hlh-30 vs WT are not converted.
Caveat: no fed control in this comparison, so the list mixes basal and starvation-dependent effects.

## Murphy et al. 2019 (starved L1)

Citation: Murphy JT, et al. Simple nutrients bypass the requirement for HLH-30 in coupling lysosomal nutrient sensing to survival. PLoS Biol 2019;17(5):e3000245. doi:10.1371/journal.pbio.3000245. PMC6516633. GEO GSE126061.

Source: `https://journals.plos.org/plosbiology/article/file?type=supplementary&id=10.1371/journal.pbio.3000245.s032` (S2 Table), saved as `raw/journal.pbio.3000245.s032.xlsx` (21 KB). S5 Table (`.s035`, 18 KB, refeeding clusters) is saved, not converted.

Design: N2 and hlh-30(tm1978) L1 larvae, fed vs starved 33 h, RNA-seq, n = 2. The paper states no explicit DE threshold.

- `murphy2019_hlh30_starved_l1_kegg_clusters.tsv`: heatmap groups with `group` and `kegg_category`.
  - B (38) and C (40), "Down in mutant starved": `direction = down`, lower in starved hlh-30 than in starved WT.
  - D (53), "Up in wt starved only": `direction = up`, starvation-induced only in WT (hlh-30-dependent).
  - E (54), "Up in mutant starved only": `direction = up_in_hlh30`.
  - Group A (commonly down in starvation) is left out.

Caveat: the table lists only cluster members that fall into 12 KEGG categories, so these are partial lists.
Webster AK et al. 2022 Cell Rep (PMC9608353) re-lists the same groups with WBGene IDs; not downloaded (47.7 MB, duplicate).

## Metcalf et al. 2024 (glial XBP-1s)

Citation: Metcalf MG, et al. Cell non-autonomous control of autophagy and metabolism by glial cells. iScience 2024;27(4):109354. doi:10.1016/j.isci.2024.109354. PMC10946330.

Source: `https://pmc-oa-opendata.s3.amazonaws.com/PMC10946330.1/mmc2.xlsx` (Table S2, 7.3 MB) saved as `raw/PMC10946330_mmc2.xlsx`. Table S3 (`PMC10946330_mmc3.xlsx`, 15 KB, 36 known HLH-30 targets up in glial XBP-1s) is saved, not converted.

Design: N2, glial XBP-1s, hlh-30(tm1978), glial XBP-1s;hlh-30(tm1978), day 2 adults, Kallisto + DESeq2. The paper's DEG threshold is padj < 0.05 and |log2FC| > 0.5.

- `metcalf2024_glial_xbp1s_hlh30_vs_xbp1s_de.tsv`: sheet "DEG", block "glial XBP-1s; hlh-30(tm1978) vs glial XBP-1s", padj < 0.05 and |log2FC| > 0.5. `up`/`down` = higher/lower in the hlh-30 background. The sign was checked with hlh-30 itself (log2FC -1.15). The block labels are not consistent in sign across the three blocks, so I used this internal check rather than the label.

Caveat: no plain hlh-30 vs N2 comparison; this is the effect of hlh-30 loss in the glial XBP-1s background.

## ENCODE / modERN HLH-30 ChIP-seq, late embryo

Citation: Kudron MM, et al. Genetics 2018;208(3):937-949. doi:10.1534/genetics.117.300657. ENCODE experiment ENCSR091KQX (original modENCODE submission GSE48752).

Source: `https://www.encodeproject.org/files/ENCFF574KYZ/@@download/ENCFF574KYZ.bed.gz`, optimal IDR thresholded peaks, ce11, saved as `raw/ENCFF574KYZ.bed.gz` (42 KB, 1,959 peaks). Strain OP433 (wgIs433, HLH-30::eGFP), late embryos, anti-GFP.

- `encode_modern_hlh30_late_embryo_chip_promoter_bound.tsv`: same rule and code as the L4 list in `../README.md` (summit from 1,000 bp upstream to 200 bp downstream of a protein-coding TSS). `direction = bound`.

ENCSR246PQU (L4, already in `../README.md`) and ENCSR091KQX are the only released ENCODE HLH-30 experiments.

## Found, but no usable gene table

- Ricaurte-Perez C, et al. 2025 Nat Commun, doi:10.1038/s41467-025-64832-x. Saved `raw/41467_2025_64832_MOESM5_ESM.xlsx` (16.4 MB, full DESeq2 tables) and `raw/41467_2025_64832_MOESM6_ESM.xlsx` (62 KB). The contrasts are intestinal dSVIP overexpression vs control within N2, hlh-30, daf-16 and the double mutant. There is no hlh-30 vs N2 contrast and the authors define no hlh-30-dependent set. Not converted.
- Shalash R, et al. 2025 Aging Cell, doi:10.1111/acel.70038. Saved `raw/ACEL-24-e70038-s006.xlsx` (106 KB). Sheet "HLH30 CHIP (without linc)" is 5,103 genes taken from the modENCODE HLH-30 ChIP (via Price et al. 2023), not new data. Not converted; use the ENCODE lists.
- Nakamura S, et al. 2016 Nat Commun, doi:10.1038/ncomms10944. hlh-30 and glp-1;hlh-30 RNA-seq exist, but the supplements hold only DAVID GO/KEGG terms. ArrayExpress E-MTAB-3686 has raw FASTQ only.
- Wong SQ, et al. 2023 Aging Cell, doi:10.1111/acel.13741, GEO GSE200383. N2, hlh-30(tm1978) and neuronal rescue at 20 and 37 C. GEO has per-sample counts only (39.5 MB, not downloaded); no gene list in the SI.
- Lin XX, et al. 2018 Nat Commun: raw reads only (see `../README.md`).
- Sohn 2025 hermaphrodites, Possik 2023 plain hlh-30 vs WT, Nonninger 2025 refeeding, Munoz-Barrera 2026 GEO counts: counts only. A list would need a new DE run on author counts; not done.
- Lim CY, et al. 2023 Autophagy (SAMS-1; abstract mentions ChIP-seq), doi:10.1080/15548627.2022.2068267: not open access in Europe PMC, the T&F supplement returned 403.
- Armendariz LA, et al. 2026 bioRxiv (P. aeruginosa liquid killing, HLH-30-dependent metabolic rewiring), Schmauck-Medina T, et al. 2026 bioRxiv (ATG-18), Moreno TM, et al. 2026 bioRxiv (hormetic heat), Mallick S, et al. 2026 bioRxiv (cholinergic HLH-30): bioRxiv returned HTTP 429 to both curl and the web fetcher, so full text and supplements were not read. Abstracts mention no hlh-30 transcriptome table, and GEO has no matching series. Not verified.
- Cornwell 2024 Geroscience (PMC11336136): no full text reachable; not checked.

Checked and no hlh-30 genome-wide data: Goswamy 2023 Front Immunol (nhr-42 only), Tsai 2021 Genetics (qPCR), Gonzalez & Irazoqui 2024, Hsieh 2024 Autophagy (PRMT-7), El-Houjeiri 2019 Cell Rep (flcn-1 and aak-1;aak-2 only), Najibi 2016 Cell Rep, Labed 2018 Immunity, DiLoreto 2025, Gupta 2019 bioRxiv (WT only), Otarigho 2026, Dennis 2026, Liu 2020 JCB, Silvestrini 2018 Cell Rep, van der Rijt 2025 (plag-15 only), Butler 2019, Lazaro-Pena 2023, Shi 2019 eLife, Chamoli 2023, Van Raamsdonk 2026, Goh 2018 (nhr-49 only), Cubillas 2026 (hizr-1 only; hlh-30 by qPCR), Zhou 2026 Cell Rep, Zhong 2025, Das 2025, Ghaddar 2023, Mony 2021, Kumsta 2017 (both papers), Murley 2025, Saha 2026, Wu 2025, Horikawa 2024, Matty 2022, Ji 2020.

## Search

Run on 2026-10-10.

- Europe PMC REST:
  - `"hlh-30" AND ("RNA-seq" OR microarray OR transcriptom* OR "differentially expressed" OR proteom*)`: 370 hits.
  - `"hlh-30" AND (RNA-seq OR microarray OR "differentially expressed")`: 232 hits.
  - `(ABSTRACT:"hlh-30" OR TITLE:"hlh-30" OR ABSTRACT:"HLH-30/TFEB")`: 153 hits; every C. elegans paper in this set was screened by title, and the likely ones were checked.
  - `"HLH-30" AND elegans AND (RNA-seq OR microarray OR transcriptome)` and `(hlh-30 OR "HLH-30/TFEB") AND (ChIP-seq OR "ChIP-chip")`: screened.
  - Citations of Lapierre 2013 (467), O'Rourke 2013 (290) and Visvikis 2014 (277), 821 distinct, intersected with the 232-hit query: 137 papers, of which 103 new ones with PMC full text were scanned for sentences with hlh-30 and a transcriptomics term.
- GEO (NCBI E-utilities, db=gds): `hlh-30[All Fields] AND "Caenorhabditis elegans"[Organism]`, `HLH-30[All Fields]`, `TFEB[All Fields] AND "Caenorhabditis elegans"[Organism]`, `tm1978[All Fields]`: 29 distinct GSE series, each checked.
- ArrayExpress/BioStudies `hlh-30`: 6 studies (E-GEOD-48751, -48752, -15762, -57739, E-MTAB-3686, E-MTAB-15698).
- ENCODE: all experiments with target hlh-30 (2 released, 2 archived).
- WormBase expression clusters: wormbase.org, rest.wormbase.org, api.wormbase.org and WormMine all returned a Cloudflare challenge to scripted requests, so expression clusters could not be queried.
  As a substitute, the WormExp v2 gene-set download (2017) was searched: its hlh-30 sets are Grove Table S8 (up 12, down 123), Grove Table S6 (PBM targets), Visvikis Table S3 and Chen 2017 Table S1, all covered here or in `../README.md`.
  Expression clusters added to WormBase after 2017 were not checked.
- bioRxiv full text: blocked (HTTP 429).
