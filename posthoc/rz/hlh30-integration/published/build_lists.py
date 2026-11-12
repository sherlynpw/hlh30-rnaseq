"""Build tidy HLH-30 gene-list TSVs from the raw supplementary files in raw/.

Run from this directory:
    uv run --with pandas --with openpyxl python build_lists.py
"""

import gzip
import re
import sqlite3
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
REPO = HERE.parents[3]
GENE_DB = REPO / "gene_db.sqlite"
GTF = next((REPO / "work").glob("*/*/Caenorhabditis_elegans.WBcel235.115.gtf"))

# Peak summit counts as promoter-bound when it lies this far up/downstream of any
# protein-coding transcript start (strand-aware).
TSS_UPSTREAM = 1000
TSS_DOWNSTREAM = 200


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
    """Return (gene_id, current gene_name, method) for a gene or sequence name."""
    key = str(raw).strip().lower()
    for method, index in (("gene_name", BY_NAME), ("sequence_name", BY_SEQ), ("other_names", BY_ALIAS)):
        if key in index:
            return (*index[key], method)
    return ("", "", "unmapped")


def with_mapping(df, name_col):
    mapped = df[name_col].map(map_name)
    df.insert(0, "gene_id", [m[0] for m in mapped])
    df["wb_gene_name"] = [m[1] for m in mapped]
    df["map_method"] = [m[2] for m in mapped]
    return df


def write(df, name):
    df.to_csv(HERE / f"{name}.tsv", sep="\t", index=False)
    n = len(df)
    n_mapped = (df["gene_id"] != "").sum()
    print(f"{name}\t{n}\t{n_mapped}")


def visvikis():
    xlsx = RAW / "1-s2.0-S1074761314001794-mmc2.xlsx"
    two = ["wt_ecoli_fpkm", "wt_saureus_fpkm", "wt_fc", "hlh30_ecoli_fpkm", "hlh30_saureus_fpkm", "hlh30_fc", "wt_pvalue", "hlh30_pvalue"]
    specs = {
        "visvikis2014_saureus_induced_wt": ("Table S2", ["wt_ecoli_fpkm", "wt_saureus_fpkm", "wt_fc", "wt_pvalue"], "up"),
        "visvikis2014_hlh30_dependent_saureus_induced": ("Table S3", two, "up"),
        "visvikis2014_hlh30_independent_saureus_induced": ("Table S5", two, "up"),
    }
    for name, (sheet, cols, direction) in specs.items():
        df = pd.read_excel(xlsx, sheet_name=sheet, header=2).dropna(subset=["Gene"])
        df = df.iloc[:, 1:]
        df.columns = ["source_gene"] + cols
        df = with_mapping(df, "source_gene")
        df.insert(2, "direction", direction)
        df = df.rename(columns={"wb_gene_name": "gene_name"})
        df = df[["gene_id", "gene_name", "source_gene", "direction", *cols, "map_method"]]
        write(df, name)

    raw = pd.read_excel(xlsx, sheet_name="Table S6", header=None)
    start = raw.index[raw[1].astype(str).str.strip() == "Gene"][0] + 1
    df = raw.loc[start:, [1]].dropna().rename(columns={1: "source_gene"})
    df = with_mapping(df, "source_gene").rename(columns={"wb_gene_name": "gene_name"})
    df["direction"] = "up"
    write(df[["gene_id", "gene_name", "source_gene", "direction", "map_method"]], "visvikis2014_hlh30_dependent_ebox")


def orourke():
    df = pd.read_excel(RAW / "41556_2013_BFncb2741_MOESM32_ESM.xlsx", header=None, skiprows=4).dropna(how="all")
    df.columns = ["source_gene", "wt_fasting_fc_mean_sem", "hlh30_fasting_fc_mean_sem"]
    df = with_mapping(df, "source_gene").rename(columns={"wb_gene_name": "gene_name"})
    df["direction"] = "up"
    write(df[["gene_id", "gene_name", "source_gene", "direction", "wt_fasting_fc_mean_sem", "hlh30_fasting_fc_mean_sem", "map_method"]], "orourke2013_hlh30_dependent_autophagy_qpcr")


def lapierre():
    # Table S1 of the SI PDF, transcribed by hand (gene, cosmid, number of CLEAR-like sites listed).
    rows = [
        ("lgg-1", "C32D5.9", 2), ("atg-18", "F41E6.13", 3), ("atg-9", "T22H9.2", 3), ("sqst-1", "T12G3.1", 2),
        ("vps-11", "R06F6.2", 0), ("vps-18", "W06B4.3", 1), ("lmp-1", "C03B1.12", 3), ("vha-15", "T14F9.1", 1),
        ("vha-16", "C30F8.2", 0), ("vha-17", "F49C12.13", 1), ("sul-1", "K09C4.8", 3), ("sul-2", "D1014.1", 1),
        ("sul-3", "C54D2.4", 5), ("ctsa-1", "C08H9.1", 2), ("cpr-1", "C52E4.1", 3), ("asp-1", "Y39B6A.20", 0),
        ("hlh-30", "W02C12.3", 1),
    ]
    df = pd.DataFrame(rows, columns=["source_gene", "sequence_name", "n_clear_sites_7of10"])
    mapped = [map_name(s) if map_name(s)[0] else map_name(g) for g, s in zip(df.source_gene, df.sequence_name)]
    df.insert(0, "gene_id", [m[0] for m in mapped])
    df.insert(1, "gene_name", [m[1] for m in mapped])
    df["direction"] = "NA"
    df["map_method"] = ["sequence_name" if map_name(s)[0] else m[2] for s, m in zip(df.sequence_name, mapped)]
    write(df[["gene_id", "gene_name", "source_gene", "sequence_name", "direction", "n_clear_sites_7of10", "map_method"]], "lapierre2013_clear_motif_candidates")


def harvald():
    xlsx = RAW / "1-s2.0-S2405471217302326-mmc2.xlsx"
    out = []
    for sheet, direction in (("Grove et al. HLH-30 targets_dec", "down"), ("Grove et al. HLH-30 targets_inc", "up")):
        df = pd.read_excel(xlsx, sheet_name=sheet, header=1)
        df = df[df["From"].astype(str).str.startswith("WBGene")]
        df = df.iloc[:, [0, 1, 4, 8, 9, 13, 14]]
        df.columns = ["gene_id", "uniprot", "source_gene", "protein_B", "protein_score", "mrna_B", "mrna_score"]
        df["direction"] = direction
        out.append(df)
    df = pd.concat(out).drop_duplicates()
    df["gene_name"] = df["gene_id"].map(lambda g: next((v[1] for v in BY_NAME.values() if v[0] == g), ""))
    df["map_method"] = "source"
    write(df[["gene_id", "gene_name", "source_gene", "uniprot", "direction", "protein_B", "protein_score", "mrna_B", "mrna_score", "map_method"]], "harvald2017_known_hlh30_targets")

    df = pd.read_excel(xlsx, sheet_name="Proteins B score > 1.85 _ CLEAR", header=None, skiprows=1)
    df = df[[1]].dropna().rename(columns={1: "source_gene"})
    df = with_mapping(df, "source_gene").rename(columns={"wb_gene_name": "gene_name"})
    df["direction"] = "NA"
    write(df[["gene_id", "gene_name", "source_gene", "direction", "map_method"]], "harvald2017_starvation_hlh30_protein_highB_clear")


def grove():
    df = pd.read_excel(RAW / "1-s2.0-S0092867409005194-mmc7.xls", header=1)
    df = df[["HLH-30"]].dropna().rename(columns={"HLH-30": "gene_id"})
    names = {v[0]: v[1] for v in BY_SEQ.values()}
    df["gene_name"] = df["gene_id"].map(names).fillna("")
    df["source_gene"] = df["gene_id"]
    df["direction"] = "NA"
    df["map_method"] = "source"
    write(df[["gene_id", "gene_name", "source_gene", "direction", "map_method"]], "grove2009_hlh30_pbm_predicted_targets")


def encode_chip():
    tss = []
    with open(GTF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "transcript" or 'gene_biotype "protein_coding"' not in f[8]:
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            gname = re.search(r'gene_name "([^"]+)"', f[8])
            pos = int(f[3]) if f[6] == "+" else int(f[4])
            tss.append(("chr" + f[0] if f[0] != "MtDNA" else "chrM", pos, f[6], gid, gname.group(1) if gname else ""))
    tss = pd.DataFrame(tss, columns=["chrom", "tss", "strand", "gene_id", "gene_name"]).drop_duplicates()

    with gzip.open(RAW / "ENCFF793PTZ.bed.gz", "rt") as fh:
        peaks = pd.read_csv(fh, sep="\t", header=None)
    peaks = peaks.iloc[:, [0, 1, 6, 8, 9]]
    peaks.columns = ["chrom", "start", "signal", "neglog10_q", "summit_offset"]
    peaks["summit"] = peaks["start"] + peaks["summit_offset"] + 1

    hits = []
    for chrom, p in peaks.groupby("chrom"):
        t = tss[tss.chrom == chrom]
        for row in p.itertuples():
            rel = row.summit - t.tss.values
            rel = rel * (t.strand.values == "+") - rel * (t.strand.values == "-")
            ok = (rel >= -TSS_UPSTREAM) & (rel <= TSS_DOWNSTREAM)
            for (gid, gname), d in zip(t[ok][["gene_id", "gene_name"]].itertuples(index=False), rel[ok]):
                hits.append((gid, gname, int(d), row.signal))
    hits = pd.DataFrame(hits, columns=["gene_id", "gene_name", "summit_to_tss", "signal"])
    hits["absd"] = hits.summit_to_tss.abs()
    best = hits.sort_values(["gene_id", "absd"]).groupby("gene_id", as_index=False).first()
    best["direction"] = "bound"
    best["source_gene"] = best["gene_name"]
    best["map_method"] = "gtf"
    best = best.sort_values("signal", ascending=False)
    write(best[["gene_id", "gene_name", "source_gene", "direction", "summit_to_tss", "signal", "map_method"]], "encode_modern_hlh30_l4_chip_promoter_bound")
    print(f"# peaks {len(peaks)}; genes {len(best)}")


if __name__ == "__main__":
    visvikis()
    orourke()
    lapierre()
    harvald()
    grove()
    encode_chip()
