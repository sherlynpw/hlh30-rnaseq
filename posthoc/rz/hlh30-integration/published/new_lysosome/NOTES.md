# Core lysosome genes (Sun et al. 2020)

Source: Sun Y, Li M, Zhao D, Li X, Yang C, Wang X. Lysosome activity is modulated by multiple longevity pathways and is important for lifespan extension in *C. elegans*. eLife. 2020;9:e55745. doi:10.7554/eLife.55745

The list is their Supplementary file 1, "The 85 lysosome-related genes analyzed by RT-PCR", in `raw/sun2020_elife-55745-supp1-v2.docx`.
It has five categories: lysosomal membrane proteins (11), V0 subunits (11), V1 subunits (9), proteases (cathepsins) (26) and non-protease hydrolases (28).
The docx has one empty row in the V1 block, which was dropped.

| file | meaning |
|---|---|
| sun2020_lysosome_genes_raw.tsv | the docx table as typed (category, gene, description) |
| sun2020_lysosome_genes.tsv | the same rows with WBGene IDs (85 rows, 85 distinct IDs) |

Names map to WBGene IDs as in `../new_autophagy/build_lists.py`: a case-insensitive match against `gene_name`, then `sequence_name`, then `other_names` in `gene_db.sqlite`.
One name needs a fixed target, listed in `ALIASES` in `build_list.py`: `ctsa-1` maps to C08H9.1 (now ctsa-3.2), because the current alias ctsa-1 points to K10B2.2, which the table also lists by cosmid.

```
uv run --with pandas python build_list.py
```

Used by fig 04d (`../../published-overlap.R --paper`).
