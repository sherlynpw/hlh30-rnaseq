"""Map the 85 lysosome genes of Sun et al. 2020 (eLife, Supplementary file 1) to WBGene IDs.

Input:  sun2020_lysosome_genes_raw.tsv (category, gene, description), typed out of
        raw/sun2020_elife-55745-supp1-v2.docx
Output: sun2020_lysosome_genes.tsv

Run from this directory:
    uv run --with pandas python build_list.py
"""

import sqlite3
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
GENE_DB = HERE.parents[4] / "gene_db.sqlite"

# Author spellings that need a fixed target before lookup.
ALIASES = {
    # The old name ctsa-1 is now an alias of ctsa-4.1 (K10B2.2), which the table
    # also lists by cosmid, so that cannot be the gene meant. Of the two CTSA
    # genes the table lacks (C08H9.1, F13D12.6), C08H9.1 carried the name ctsa-1
    # (Lapierre 2013 Table S1; now ctsa-3.2).
    "ctsa-1": "C08H9.1",
}


def load_name_index():
    con = sqlite3.connect(GENE_DB)
    anno = pd.read_sql(
        "select wormbase_gene, gene_name, sequence_name, other_names, status from annotations",
        con,
    )
    con.close()
    anno = anno.sort_values("status", key=lambda s: s != "Live")
    by_name, by_seq, by_alias = {}, {}, {}
    for row in anno.itertuples():
        if row.gene_name:
            by_name.setdefault(row.gene_name.lower(), (row.wormbase_gene, row.gene_name))
        if row.sequence_name:
            by_seq.setdefault(row.sequence_name.lower(), (row.wormbase_gene, row.gene_name))
        for alias in (row.other_names or "").split(";"):
            alias = alias.strip().lower()
            if alias:
                by_alias.setdefault(alias, (row.wormbase_gene, row.gene_name))
    return by_name, by_seq, by_alias


BY_NAME, BY_SEQ, BY_ALIAS = load_name_index()


def map_name(raw):
    key = str(raw).strip().lower()
    aliased = key in ALIASES
    key = ALIASES.get(key, key).lower()
    for method, index in (("gene_name", BY_NAME), ("sequence_name", BY_SEQ), ("other_names", BY_ALIAS)):
        if key in index:
            return (*index[key], "alias" if aliased else method)
    return ("", "", "unmapped")


raw = pd.read_csv(HERE / "sun2020_lysosome_genes_raw.tsv", sep="\t", keep_default_na=False)
mapped = raw["gene"].map(map_name)
out = pd.DataFrame({
    "gene_id": [m[0] for m in mapped],
    "gene_name": [m[1] for m in mapped],
    "source_gene": raw["gene"],
    "category": raw["category"],
    "description": raw["description"],
    "map_method": [m[2] for m in mapped],
})
out.to_csv(HERE / "sun2020_lysosome_genes.tsv", sep="\t", index=False)
print(f"{len(out)} rows, {(out.gene_id != '').sum()} mapped, {out.loc[out.gene_id != '', 'gene_id'].nunique()} distinct IDs")
print(out[out.map_method != "gene_name"].to_string(index=False))
