# hlh-30 RNA-seq analysis

Analysis code for the *hlh-30(tm1978)* RNA-seq study of *C. elegans* at L1, L4, Day 3 and Day 6 of adulthood.
The RNA-seq data are in NCBI GEO under accession GSE354145.

## Layout

- `main.nf`, `modules/`, `subworkflows/`, `bin/` and `conf/` are the Nextflow pipeline.
  It counts reads with featureCounts, removes unwanted variation with RUVSeq, tests differential expression with DESeq2, clusters genes and runs the enrichment tests.
- `samples/` holds the sample sheets and the run configurations of the paper.
- `posthoc/rz/` holds the scripts that draw the figure panels and build the supplementary files from the pipeline output.
  `posthoc/rz/paper-figures/PROVENANCE.md` maps each panel to its script and input.

## Pipeline runs

| Run | Samples | Used for |
|---|---|---|
| `VZRJXQ-ruv-rbe` | WT and *hlh-30* at L4, Day 3 and Day 6 | Fig. 1B, 1C, 1E, 1F, 3A to 3D, 4, S1A, S1B and S2A to S2E |
| `VZRJXQ_no_l4-ruv-rbe` | WT and *hlh-30* at Day 3 and Day 6 | Fig. 1D, 3E to 3H and S1C |
| `hlh30-all-ruv` | WT and *hlh-30* at L1, L4, Day 3 and Day 6 | Fig. 9, S2F, S3, S5C and S7 to S11 |
| `pca-hlh30-all-ruv` | as `hlh30-all-ruv` | Fig. 9A and 9C |
| `pca-hlh30-l1` | WT and *hlh-30* at L1 | Fig. 9B |

Each run reads Plasmidsaurus deduplicated BAM files, so the alignment step is skipped:

```bash
./run.sh --params samples/VZRJXQ-ruv-rbe.yaml --profile aws --skip-alignment
```

## Figures

Each folder in `posthoc/rz/` has a Makefile that downloads its inputs from the pipeline output and renders its panels.
`make -C posthoc/rz/paper-figures` collects the panels under their figure names.
