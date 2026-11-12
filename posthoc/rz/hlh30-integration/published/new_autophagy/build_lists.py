"""Build tidy autophagy gene-list TSVs for the HLH-30 integration.

Inputs:
    evidence_curated.tsv            gene-level claims typed in from the papers (see NOTES.md)
    raw/europepmc_xml/PMC12005284.xml  Zhang & Melendez 2025 Genetics, Table 1
    raw/europepmc_xml/PMC4502811.xml   Zhang et al. 2015 Autophagy, Table 1
    raw/go0006914_orgCeegdb.tsv     from export_go.R

Run from this directory:
    Rscript export_go.R
    uv run --with pandas python build_lists.py
"""

import re
import sqlite3
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
REPO = HERE.parents[4]
GENE_DB = REPO / "gene_db.sqlite"

# Author spellings that need a fixed target before lookup.
ALIASES = {
    # Lapierre 2013 figure label; their Table S1 gives cosmid C08H9.1 (now ctsa-3.2).
    # The old name ctsa-1 is now an alias of ctsa-4.1, so map by cosmid. Chamoli 2023
    # writes "ctsa" with no cosmid and is left unmapped.
    "ctsa*": "C08H9.1",
    "y34ba.2": "Y34B4A.2",  # Zhang & Melendez 2025 Table 1 typo for the UVRAG homolog
    "impk-1": "ipmk-1",  # Table 1 "impk-1/IMPK"; the worm IPMK homolog is ipmk-1 (ZK795.1)
    "fundc-1": "fndc-1",  # Table 1 "fundc-1/FUNDC1A"; the worm FUNDC1 homolog is fndc-1 (T06D8.7)
    "rub-1": "Y56A3A.16",  # Table 1 "rub-1/RUBCN"; WormBase names Y56A3A.16 as the RUBCN ortholog
    "allo-1a": "allo-1",  # isoforms of allo-1 (R102.5)
    "allo-1b": "allo-1",
}


def load_name_index():
    con = sqlite3.connect(GENE_DB)
    anno = pd.read_sql(
        "select wormbase_gene, gene_name, sequence_name, other_names, status from annotations",
        con,
    )
    con.close()
    anno = anno.sort_values("status", key=lambda s: s != "Live")
    by_name, by_seq, by_alias, by_id = {}, {}, {}, {}
    for row in anno.itertuples():
        by_id.setdefault(row.wormbase_gene, row.gene_name)
        if row.gene_name:
            by_name.setdefault(row.gene_name.lower(), (row.wormbase_gene, row.gene_name))
        if row.sequence_name:
            by_seq.setdefault(row.sequence_name.lower(), (row.wormbase_gene, row.gene_name))
        for alias in (row.other_names or "").split(";"):
            alias = alias.strip().lower()
            if alias:
                by_alias.setdefault(alias, (row.wormbase_gene, row.gene_name))
    return by_name, by_seq, by_alias, by_id


BY_NAME, BY_SEQ, BY_ALIAS, BY_ID = load_name_index()


def map_name(raw):
    """Return (gene_id, current gene_name, method) for a gene or sequence name."""
    key = str(raw).strip().lower()
    aliased = key in ALIASES
    key = ALIASES.get(key, key).lower()
    for method, index in (("gene_name", BY_NAME), ("sequence_name", BY_SEQ), ("other_names", BY_ALIAS)):
        if key in index:
            return (*index[key], "alias" if aliased else method)
    return ("", "", "unmapped")


def with_mapping(df, name_col="source_gene"):
    mapped = df[name_col].map(map_name)
    df.insert(0, "gene_id", [m[0] for m in mapped])
    df.insert(1, "gene_name", [m[1] for m in mapped])
    df["map_method"] = [m[2] for m in mapped]
    return df


def write(df, name):
    df.to_csv(HERE / f"{name}.tsv", sep="\t", index=False)
    print(f"{name}\t{len(df)}\t{(df['gene_id'] != '').sum()}\t{df.loc[df.gene_id != '', 'gene_id'].nunique()}")


def text(e):
    return re.sub(r"\s+", " ", "".join(e.itertext())).strip() if e is not None else ""


def table_rows(pmcid):
    root = ET.parse(RAW / "europepmc_xml" / f"{pmcid}.xml").getroot()
    tw = next(root.iter("table-wrap"))
    return [[text(c) for c in tr if c.tag in ("td", "th")] for tr in tw.iter("tr")]


def split_genes(cell):
    cell = cell.replace("/TOR", "")
    return [g.strip() for g in re.split(r",|/(?=[a-z]{3,4}-)|\band\b", cell) if g.strip()]


CORE_CATEGORIES = {
    "Atg1/ULK complex": "core_machinery",
    "Class III (PI(3)K) Vps34 complex": "core_machinery",
    "PI3P binding complex": "core_machinery",
    "Atg8 ubiquitin-like conjugation system": "core_machinery",
    "Atg12 ubiquitin-like conjugation system": "core_machinery",
    "Autophagosome maturation": "maturation_fusion",
    "Autophagy-related": "autophagy_related_or_regulator",
    "Regulate autophagy gene expression": "transcriptional_regulator",
    "Selective autophagy receptor for ubiquitinated cargo": "selective_autophagy",
    "Selective autophagy– mitophagy": "selective_autophagy",
    "Selective autophagy– allophagy": "selective_autophagy",
    "HOPS complex, for vesicle fusion": "maturation_fusion",
}


def zhang_melendez_2025():
    rows = table_rows("PMC12005284")[1:]
    out, category = [], ""
    for r in rows:
        if len(r) == 4:
            category = r[0] or category
            r = r[1:]
        yeast, worm, mammal = r
        for g in split_genes(worm):
            out.append((g, category, yeast, mammal))
    df = pd.DataFrame(out, columns=["source_gene", "category", "yeast_homolog", "mammalian_homolog"])
    df["role_class"] = df["category"].map(CORE_CATEGORIES)
    z15 = zhang_2015_genes()
    df = with_mapping(df)
    df["in_zhang2015_table1"] = df["gene_id"].isin(z15)
    df = df.drop_duplicates(subset=["gene_id", "source_gene", "category"])
    write(df[["gene_id", "gene_name", "source_gene", "category", "role_class", "yeast_homolog", "mammalian_homolog", "in_zhang2015_table1", "map_method"]], "core_autophagy_zhang_melendez2025")
    return df


def zhang_2015_genes():
    ids = set()
    for r in table_rows("PMC4502811")[2:]:
        g = r[2] if len(r) >= 3 else ""
        if g:
            ids.add(map_name(g)[0])
    ids.discard("")
    return ids


def go_autophagy(core):
    df = pd.read_csv(RAW / "go0006914_orgCeegdb.tsv", sep="\t", comment="#")
    df["gene_name"] = df["gene_id"].map(BY_ID).fillna("")
    df["source_gene"] = df["symbol"]
    df["in_core_list"] = df["gene_id"].isin(core["gene_id"])
    df["map_method"] = "source"
    write(df[["gene_id", "gene_name", "source_gene", "evidence_codes", "iea_only", "in_core_list", "map_method"]], "go0006914_autophagy_orgCeegdb")
    return df


EVIDENCE_ORDER = ["expr_reduced_hlh30", "induced_hlh30_oe", "chip_binding", "motif_only", "not_dependent"]


def evidence(core, go):
    df = pd.read_csv(HERE / "evidence_curated.tsv", sep="\t", dtype=str).fillna("")
    df = with_mapping(df)
    df["in_core_list"] = df["gene_id"].isin(core["gene_id"])
    df["in_go0006914"] = df["gene_id"].isin(go["gene_id"])
    role = core.drop_duplicates("gene_id").set_index("gene_id")["role_class"]
    df["gene_class"] = df["gene_id"].map(role).fillna("lysosome_or_other")
    cols = ["gene_id", "gene_name", "source_gene", "gene_class", "evidence", "hlh30_perturbation", "condition", "stage_tissue", "assay", "result", "paper", "figure", "pmcid", "doi", "in_core_list", "in_go0006914", "quote", "notes", "map_method"]
    df = df[cols]
    write(df, "hlh30_autophagy_evidence")

    pos = df[df.evidence != "not_dependent"]
    rows = []
    for gid, d in pos.groupby("gene_id"):
        row = {"gene_id": gid, "gene_name": d.gene_name.iloc[0], "source_gene": ";".join(sorted(set(d.source_gene))), "gene_class": d.gene_class.iloc[0]}
        for ev in EVIDENCE_ORDER[:-1]:
            papers = sorted(set(d.loc[d.evidence == ev, "paper"]))
            row[f"n_papers_{ev}"] = len(papers)
            row[f"papers_{ev}"] = "; ".join(papers)
        neg = sorted(set(df.loc[(df.gene_id == gid) & (df.evidence == "not_dependent"), "paper"]))
        row["papers_not_dependent"] = "; ".join(neg)
        row["strongest_evidence"] = next(ev for ev in ["chip_binding", "expr_reduced_hlh30", "induced_hlh30_oe", "motif_only"] if row[f"n_papers_{ev}"])
        row["n_papers_functional"] = len(set(d.loc[d.evidence != "motif_only", "paper"]))
        row["in_core_list"] = d.in_core_list.iloc[0]
        row["in_go0006914"] = d.in_go0006914.iloc[0]
        row["map_method"] = d.map_method.iloc[0]
        rows.append(row)
    summary = pd.DataFrame(rows).sort_values(["n_papers_functional", "gene_name"], ascending=[False, True])
    write(summary, "hlh30_autophagy_genes_summary")


if __name__ == "__main__":
    core = zhang_melendez_2025()
    go = go_autophagy(core)
    evidence(core, go)
