"""Build tidy TSVs for the additional HLH-30 gene lists in raw/.

Same conventions as ../build_lists.py. Run from this directory:
    uv run --with pandas --with openpyxl --with xlrd python build_new.py
"""

import gzip
import re
import sqlite3
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
REPO = HERE.parents[4]
GENE_DB = REPO / "gene_db.sqlite"
GTF = next((REPO / "work").glob("*/*/Caenorhabditis_elegans.WBcel235.115.gtf"))

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
    by_id = dict(zip(anno.wormbase_gene, anno.gene_name))
    return by_name, by_seq, by_alias, by_id


BY_NAME, BY_SEQ, BY_ALIAS, BY_ID = load_name_index()


def map_name(raw):
    """Return (gene_id, current gene_name, method) for a gene or sequence name."""
    key = str(raw).strip().lower()
    for method, index in (("gene_name", BY_NAME), ("sequence_name", BY_SEQ), ("other_names", BY_ALIAS)):
        if key in index:
            return (*index[key], method)
    return ("", "", "unmapped")


def map_column(df, name_col):
    """Add gene_id, gene_name and map_method from a gene/sequence name column."""
    mapped = df[name_col].map(map_name)
    df["gene_id"] = [m[0] for m in mapped]
    df["gene_name"] = [m[1] for m in mapped]
    df["map_method"] = [m[2] for m in mapped]
    return df


def map_ids(df, id_col, fallback_col=None):
    """Use author WBGene IDs; fall back to a name column for rows without one."""
    ids = df[id_col].astype(str).str.strip()
    df["gene_id"] = ids.where(ids.str.startswith("WBGene"), "")
    df["gene_name"] = df["gene_id"].map(BY_ID).fillna("")
    df["map_method"] = df["gene_id"].map(lambda g: "source" if g else "unmapped")
    if fallback_col is not None:
        miss = df["gene_id"] == ""
        mapped = df.loc[miss, fallback_col].map(map_name)
        df.loc[miss, "gene_id"] = [m[0] for m in mapped]
        df.loc[miss, "gene_name"] = [m[1] for m in mapped]
        df.loc[miss, "map_method"] = [m[2] for m in mapped]
    return df


def write(df, name, extra):
    cols = ["gene_id", "gene_name", "source_gene", "direction", *extra, "map_method"]
    df = df[cols]
    df.to_csv(HERE / f"{name}.tsv", sep="\t", index=False)
    n_mapped = (df["gene_id"] != "").sum()
    n_ids = df.loc[df["gene_id"] != "", "gene_id"].nunique()
    by_dir = ", ".join(f"{k}={v}" for k, v in df["direction"].value_counts().items())
    print(f"{name}\t{len(df)}\t{n_mapped}\t{n_ids}\t{by_dir}")


def grove_microarray():
    raw = pd.read_excel(RAW / "1-s2.0-S0092867409005194-mmc9.xls", header=None)
    first = raw[0].astype(str)
    up_start = first.index[first.str.startswith("Genes with increased")][0]
    down_start = first.index[first.str.startswith("Genes with decreased")][0]
    rows = []
    for start, stop, direction in ((up_start, down_start, "up"), (down_start, len(raw), "down")):
        block = raw.iloc[start + 1 : stop]
        block = block[block[0].astype(str).str.startswith("WBGene")].copy()
        block["direction"] = direction
        rows.append(block)
    df = pd.concat(rows)
    df.columns = ["source_id", "source_gene", "sequence_name", "logFC", "pvalue", "padj", "direction"]
    df = map_ids(df, "source_id", fallback_col="sequence_name")
    write(df, "grove2009_hlh30_microarray_de", ["sequence_name", "logFC", "pvalue", "padj"])


def sohn_males():
    df = pd.read_excel(RAW / "kaup_a_2375779_sm6158.xlsx", sheet_name="Fig. 5A", header=3).dropna(subset=["Gene ID"])
    df.columns = [
        "source_gene", "source_id",
        "male_vs_herm_log2fc", "male_vs_herm_lfcse", "male_vs_herm_pvalue", "male_vs_herm_padj",
        "hlh30_male_vs_wt_male_log2fc", "hlh30_male_vs_wt_male_lfcse", "hlh30_male_vs_wt_male_pvalue", "hlh30_male_vs_wt_male_padj",
    ]
    df["direction"] = "down"
    df = map_ids(df, "source_id", fallback_col="source_gene")
    extra = ["hlh30_male_vs_wt_male_log2fc", "hlh30_male_vs_wt_male_padj", "male_vs_herm_log2fc", "male_vs_herm_padj"]
    write(df, "sohn2025_hlh30_down_in_males", extra)


def chen_cry5b():
    df = pd.read_excel(RAW / "2016AUTO0201R2-s03.xlsx", header=1).dropna(subset=["Wormbase ID"])
    df.columns = ["source_id", "source_gene", "sequence_name", "log2fc_normalized", "pvalue"]
    df = map_ids(df, "source_id", fallback_col="sequence_name")
    df["direction"] = "up"
    df = df.sort_values("pvalue").drop_duplicates("source_id")
    write(df, "chen2017_hlh30_dependent_cry5b_induced", ["sequence_name", "log2fc_normalized", "pvalue"])


def wani_saureus():
    xlsx = RAW / "elife-62775-supp3-v3.xlsx"
    ov = pd.read_excel(xlsx, sheet_name="Overlap_UP_WT_v_hlh-30(-)")
    wt_de = pd.read_excel(RAW / "elife-62775-supp1-v3.xlsx", sheet_name="up+down_WT_Starv_v_Saur", header=1)
    wt_de = wt_de.dropna(subset=["ID"]).set_index("Gene")
    specs = {
        "wani2021_hlh30_dependent_saureus_induced": "Wild type ONLY",
        "wani2021_hlh30_independent_saureus_induced": "Wild type & hlh-30(-)",
    }
    for name, col in specs.items():
        df = ov[[col]].dropna().rename(columns={col: "source_gene"})
        df["source_gene"] = df["source_gene"].astype(str).str.strip()
        df["source_id"] = df["source_gene"].map(wt_de["ID"]).fillna("")
        df["wt_log2fc"] = df["source_gene"].map(wt_de["log2FoldChange"])
        df["wt_padj"] = df["source_gene"].map(wt_de["padj"])
        df = map_ids(df, "source_id", fallback_col="source_gene")
        df["direction"] = "up"
        write(df, name, ["wt_log2fc", "wt_padj"])


def hsieh_r38():
    import docx

    table = docx.Document(RAW / "kaup_a_2735736_sm2037.docx").tables[0]
    rows = []
    for r in table.rows[1:]:
        cells = [c.text.strip() for c in r.cells]
        for i in (0, 3):
            if cells[i + 1].startswith("WBGene"):
                rows.append((int(cells[i]), cells[i + 1], cells[i + 2]))
    df = pd.DataFrame(rows, columns=["rank", "source_id", "source_gene"])
    df = df.drop_duplicates().sort_values("rank")
    df = map_ids(df, "source_id", fallback_col="source_gene")
    df["direction"] = "up"
    write(df, "hsieh2026_hlh30_r38_dependent_cry5b_induced", ["rank"])


def sign_direction(lfc):
    return lfc.map(lambda x: "up" if x > 0 else "down")


def gerisch_ard():
    df = pd.read_csv(RAW / "GSE126632_gene.expression.tsv.gz", sep="\t")
    df = df.rename(columns={"log2(fold_change)": "log2fc", "p_value": "pvalue", "q_value": "qvalue"})
    for cond in ("ARD", "AL"):
        d = df[(df.sample_1 == f"N2_{cond}") & (df.sample_2 == f"h30_{cond}") & (df.significant == "yes")].copy()
        d["source_gene"] = d["gene"]
        d["direction"] = sign_direction(d["log2fc"])
        d["paper_go_cutoff"] = (d["log2fc"].abs() >= 2) & (d["pvalue"] < 0.05)
        d = map_ids(d, "ensembl_gene_id", fallback_col="gene").sort_values("qvalue")
        write(d, f"gerisch2020_hlh30_vs_wt_{cond.lower()}_de", ["log2fc", "pvalue", "qvalue", "paper_go_cutoff"])


def possik_pgph2():
    df = pd.read_excel(RAW / "41467_2023_40857_MOESM11_ESM.xlsx", header=2).dropna(subset=["Gene ID"])
    df = df.rename(columns={"Gene ID": "sequence_name", "Gene name": "source_gene"})
    mapped = [map_name(s) if map_name(s)[0] else map_name(g) for s, g in zip(df.sequence_name, df.source_gene)]
    df["gene_id"] = [m[0] for m in mapped]
    df["gene_name"] = [m[1] for m in mapped]
    df["map_method"] = [m[2] for m in mapped]
    df["direction"] = "up"
    write(df, "possik2023_hlh30_dependent_pgph2_induced", ["sequence_name"])


def shioda_neuronal_rnai():
    df = pd.read_excel(RAW / "GSE234429_gTU_luc_vs_gTU_hlh-30i.xlsx")
    df = df[df["FDR"] < 0.05].rename(columns={"Geneid": "source_gene", "logFC": "log2fc", "p.value": "pvalue", "FDR": "fdr"})
    df["direction"] = sign_direction(df["log2fc"])
    df = map_column(df, "source_gene")
    write(df, "shioda2023_glp1_neuronal_hlh30_rnai_de", ["log2fc", "pvalue", "fdr"])


def lin2024_abeta():
    df = pd.read_csv(RAW / "GSE272146_GMC101mutvsGMC101_deg.txt.gz", sep="\t")
    df = df[(df["padj"] < 0.05) & (df["log2FoldChange"].abs() >= 1)].copy()
    df = df.rename(columns={"gene_id": "source_gene", "log2FoldChange": "log2fc"})
    df["direction"] = sign_direction(df["log2fc"])
    df = map_column(df, "source_gene")
    write(df, "lin2024_gmc101_hlh30_ko_de", ["log2fc", "pvalue", "padj"])


def nonninger_ard():
    df = pd.read_excel(RAW / "43587_2025_911_MOESM12_ESM.xlsx", sheet_name="Extended Figure 5a").dropna(subset=["Gene"])
    df.columns = ["source_gene", "ard_log2fc_hlh30_vs_n2", "ard_log2fc_hlh30daf1_vs_hlh30"]
    df["direction"] = "NA"
    df = map_column(df, "source_gene").sort_values("ard_log2fc_hlh30_vs_n2")
    write(df, "nonninger2025_hlh30_vs_wt_ard_log2fc_all_genes", ["ard_log2fc_hlh30_vs_n2", "ard_log2fc_hlh30daf1_vs_hlh30"])


def munoz_barrera_l1_starved():
    df = pd.read_excel(RAW / "gkag949_File_S2_RNAseq_Rev.xlsx", sheet_name="deseq2_hlh-30_wt")
    df = df.rename(columns={df.columns[0]: "source_id", "log2FoldChange": "log2fc", "baseMean": "base_mean"})
    df = df[df["padj"] < 0.05].copy()
    df["source_gene"] = df["source_id"]
    df["direction"] = sign_direction(df["log2fc"])
    df = map_ids(df, "source_id").sort_values("padj")
    write(df, "munozbarrera2026_hlh30_vs_wt_starved_l1_de", ["base_mean", "log2fc", "pvalue", "padj"])


def murphy_starved_l1():
    raw = pd.read_excel(RAW / "journal.pbio.3000245.s032.xlsx", header=None)
    groups = {6: ("C", "down"), 3: ("B", "down"), 9: ("D", "up"), 12: ("E", "up_in_hlh30")}
    rows = []
    for col, (group, direction) in groups.items():
        block = raw.iloc[5:, [col, col + 1]].dropna(subset=[col])
        for gene, category in block.itertuples(index=False):
            rows.append((str(gene).strip(), category, group, direction))
    df = pd.DataFrame(rows, columns=["source_gene", "kegg_category", "group", "direction"])
    df = df.sort_values("group").drop_duplicates(["source_gene", "group"])
    df = map_column(df, "source_gene")
    write(df, "murphy2019_hlh30_starved_l1_kegg_clusters", ["group", "kegg_category"])


def metcalf_glial_xbp1s():
    raw = pd.read_excel(RAW / "PMC10946330_mmc2.xlsx", sheet_name="DEG", header=None, skiprows=2)
    df = raw.iloc[:, [0, 1, 8, 9, 12, 13]].copy()
    df.columns = ["source_gene", "source_id", "base_mean", "log2fc", "pvalue", "padj"]
    df = df[(df["padj"] < 0.05) & (df["log2fc"].abs() > 0.5)].copy()
    df["direction"] = sign_direction(df["log2fc"])
    df = map_ids(df, "source_id", fallback_col="source_gene").sort_values("padj")
    write(df, "metcalf2024_glial_xbp1s_hlh30_vs_xbp1s_de", ["base_mean", "log2fc", "pvalue", "padj"])


def load_tss():
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
    return pd.DataFrame(tss, columns=["chrom", "tss", "strand", "gene_id", "gene_name"]).drop_duplicates()


def chip_promoter_bound(bed, name):
    tss = load_tss()
    with gzip.open(RAW / bed, "rt") as fh:
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
    write(best, name, ["summit_to_tss", "signal"])
    print(f"# {bed}: peaks {len(peaks)}; genes {len(best)}")


if __name__ == "__main__":
    print("name\trows\tmapped\tdistinct_ids\tdirections")
    grove_microarray()
    sohn_males()
    chen_cry5b()
    wani_saureus()
    hsieh_r38()
    gerisch_ard()
    possik_pgph2()
    shioda_neuronal_rnai()
    lin2024_abeta()
    nonninger_ard()
    munoz_barrera_l1_starved()
    murphy_starved_l1()
    metcalf_glial_xbp1s()
    chip_promoter_bound("ENCFF574KYZ.bed.gz", "encode_modern_hlh30_late_embryo_chip_promoter_bound")
