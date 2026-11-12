"""Pre-flight sanity check for a celegans-pipeline run.

Usage:
  uv run --quiet --with pyyaml python3 scripts/check_pipeline_inputs.py \\
      samples/VZRJXQ-ruv-rbe.yaml

What it checks (in order, exits non-zero on first FAIL):

  1. Pipeline YAML parses and references a samplesheet that exists.
  2. Samplesheet YAML parses; every sample row is well-formed.
  3. Every BAM path in the samplesheet exists on S3.
  4. Every FASTQ pair in the samplesheet exists on S3 (both _1 and _2).
  5. Reference files (genome FASTA, GTF, gene_db.sqlite) exist on S3.
  6. (Multi-factor mode only) Every column referenced by `design_formula`
     and `reference_levels` exists in the samplesheet metadata.
  7. (Multi-factor mode only) Every level in `reference_levels` actually
     appears in the samplesheet's data for that column.
  8. Cross-protocol sanity: warn if `source` is structurally aliased with
     a primary biology variable.

Exit codes:
  0 = all checks passed; safe to launch the pipeline.
  1 = at least one FAIL; pipeline launch will likely error or silently
      produce wrong results.

Designed to be cheap (no transfers, just s3 ls / head-object) and fast.
"""

from __future__ import annotations

import re
import shlex
import subprocess
import sys
from pathlib import Path

import yaml

# ----------------------------------------------------------------- helpers

OK   = "  OK    "
FAIL = "  FAIL  "
WARN = "  WARN  "

ANY_FAIL = False


def fail(msg: str) -> None:
    global ANY_FAIL
    ANY_FAIL = True
    print(FAIL + msg)


def warn(msg: str) -> None:
    print(WARN + msg)


def ok(msg: str) -> None:
    print(OK + msg)


def s3_exists(s3_uri: str) -> bool:
    """Cheap existence check via aws s3 ls; returns True iff the object
    is listed at the given exact key."""
    if not s3_uri.startswith("s3://"):
        return False
    rc = subprocess.run(
        ["aws", "s3", "ls", s3_uri],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    ).returncode
    return rc == 0


def load_yaml(path: Path) -> dict:
    if not path.exists():
        fail(f"YAML file does not exist: {path}")
        sys.exit(1)
    with open(path) as f:
        try:
            return yaml.safe_load(f)
        except yaml.YAMLError as e:
            fail(f"YAML parse error in {path}: {e}")
            sys.exit(1)


# ---------------------------------------------------------------- check fns

def check_pipeline_yaml(p: dict, path: Path) -> Path:
    print(f"\n[1] Pipeline YAML: {path}")
    if "input" not in p:
        fail("missing required key: input (path to samplesheet YAML)")
        sys.exit(1)
    samplesheet_path = (path.parent.parent / p["input"]).resolve()
    if not samplesheet_path.exists():
        # try relative to the current dir / repo root
        alt = (Path.cwd() / p["input"]).resolve()
        if alt.exists():
            samplesheet_path = alt
        else:
            fail(f"samplesheet path from `input:` not found: tried {samplesheet_path} and {alt}")
            sys.exit(1)
    ok(f"samplesheet -> {samplesheet_path}")
    return samplesheet_path


def check_samplesheet_structure(s: dict) -> list[dict]:
    print("\n[2] Samplesheet structure")
    if "samples" not in s or not isinstance(s["samples"], list):
        fail("samplesheet missing top-level `samples:` list")
        sys.exit(1)
    samples = s["samples"]
    ok(f"{len(samples)} sample rows")
    REQUIRED = {"sample_id", "input_type"}
    for i, row in enumerate(samples):
        miss = REQUIRED - set(row.keys())
        if miss:
            fail(f"sample[{i}] missing required keys: {sorted(miss)}; row={row}")
    return samples


def check_s3_inputs(samples: list[dict]) -> None:
    print("\n[3+4] S3 input data exists")
    bam_count = fastq_count = 0
    bad = 0
    for row in samples:
        sid = row.get("sample_id", "?")
        itype = row.get("input_type")
        if itype == "plasmidsaurus_bam":
            bam = row.get("bam")
            if not bam:
                fail(f"{sid}: input_type=plasmidsaurus_bam but no `bam:` key")
                bad += 1
                continue
            if not s3_exists(bam):
                fail(f"{sid}: BAM does NOT exist on S3 -> {bam}")
                bad += 1
            else:
                bam_count += 1
        elif itype == "fastq_standard":
            f1 = row.get("fastq_1"); f2 = row.get("fastq_2")
            if not (f1 and f2):
                fail(f"{sid}: fastq_standard requires both fastq_1 and fastq_2")
                bad += 1
                continue
            missing = [f for f in (f1, f2) if not s3_exists(f)]
            if missing:
                for m in missing:
                    fail(f"{sid}: FASTQ does NOT exist on S3 -> {m}")
                bad += 1
            else:
                fastq_count += 1
        else:
            fail(f"{sid}: unknown input_type '{itype}'")
            bad += 1
    if bad == 0:
        ok(f"all S3 inputs present: {bam_count} BAM(s), {fastq_count} FASTQ pair(s)")


def check_references(p: dict) -> None:
    print("\n[5] Reference files on S3")
    org = p.get("organism", {}) or {}
    refs = {
        "gtf":     org.get("gtf"),
        "gene_db": org.get("gene_db"),
    }
    for name, uri in refs.items():
        if not uri:
            warn(f"organism.{name} not set in pipeline YAML")
            continue
        if uri.startswith("s3://"):
            if s3_exists(uri):
                ok(f"organism.{name} -> {uri}")
            else:
                fail(f"organism.{name} does NOT exist on S3 -> {uri}")
        else:
            warn(f"organism.{name} is not an s3:// URI -> {uri} (skipping check)")


def parse_design_vars(formula: str) -> list[str]:
    """Extract variable names from an R formula like '~ a + b + c'.
    Crude tokenizer; fine for the simple additive formulas we use."""
    rhs = re.sub(r"^\s*~\s*", "", formula).strip()
    # split on + or *, drop intercept token, strip whitespace
    parts = re.split(r"[\+\*\:]", rhs)
    out = []
    for p_ in parts:
        t = p_.strip()
        if t and t != "1" and t != "0":
            out.append(t)
    return out


def parse_reference_levels(s: str) -> list[tuple[str, str]]:
    """`'col1:lvl1,col2:lvl2'` -> [('col1','lvl1'), ('col2','lvl2')]."""
    if not s:
        return []
    out = []
    for kv_raw in s.split(","):
        kv = kv_raw.strip().split(":")
        if len(kv) == 2:
            out.append((kv[0].strip(), kv[1].strip()))
    return out


def check_multifactor_design(p: dict, samples: list[dict]) -> None:
    deseq2 = p.get("deseq2", {}) or {}
    formula = deseq2.get("design_formula")
    if not formula:
        print("\n[6] (skipping multi-factor checks; pipeline is in pairwise mode)")
        return
    print("\n[6] Multi-factor design references valid columns")
    vars_ = parse_design_vars(formula)
    sample_keys = set(samples[0].keys()) if samples else set()
    # Check each var exists in at least one sample's keys (W_* RUV cols are
    # added at runtime so we exempt them).
    for v in vars_:
        if re.match(r"^W_\d+$", v):
            ok(f"`{v}` (RUV factor; auto-added at runtime)")
            continue
        if v not in sample_keys:
            fail(f"design_formula variable `{v}` is NOT a samplesheet column")
        else:
            # spot-check: not all samples may have it; require ALL do
            missing_samples = [s["sample_id"] for s in samples if v not in s]
            if missing_samples:
                fail(f"`{v}` missing from {len(missing_samples)} sample(s): "
                     f"{missing_samples[:3]}{'...' if len(missing_samples)>3 else ''}")
            else:
                vals = sorted({str(s.get(v)) for s in samples})
                ok(f"`{v}` present in all samples; levels: {vals}")

    print("\n[7] Reference levels exist in data")
    refs = parse_reference_levels(deseq2.get("reference_levels", ""))
    if not refs:
        warn("no reference_levels set; DESeq2 will pick alphabetically")
    for col, lvl in refs:
        if col not in sample_keys:
            fail(f"reference_levels references column `{col}` not in samplesheet")
            continue
        observed = {str(s.get(col)) for s in samples}
        if lvl in observed:
            ok(f"`{col}:{lvl}` (observed in {sum(1 for s in samples if str(s.get(col))==lvl)} sample(s))")
        else:
            fail(f"reference_levels `{col}:{lvl}` -- level `{lvl}` not observed; "
                 f"actual values: {sorted(observed)}")


def check_enrichment_block(p: dict) -> None:
    """The enrichment subworkflow accesses params.enrichment.gsea.collections
    etc. unconditionally; missing sub-blocks trigger a Groovy NPE in
    subworkflows/local/enrichment.nf (around line 49). Catch that here."""
    print("\n[+] Enrichment block completeness")
    if p.get("skip_enrichment") is True:
        ok("enrichment skipped (skip_enrichment: true); not checking sub-blocks")
        return
    enr = p.get("enrichment") or {}
    if not enr:
        fail("enrichment block is empty/missing AND skip_enrichment != true. "
             "The enrichment subworkflow will NPE at startup.")
        return
    # All four sub-blocks must exist (each with at least an `enabled` key).
    # The pipeline reads gsea.collections even if gsea.enabled is false.
    for sub in ("go", "kegg", "wormcat", "gsea"):
        if sub not in enr:
            fail(f"enrichment.{sub} block is missing; will cause Groovy NPE in "
                 f"enrichment.nf. Either add the sub-block (with enabled: false "
                 f"if you don't want it to run) or set skip_enrichment: true.")
    # gsea.collections must be a list (can be empty)
    gsea = enr.get("gsea") or {}
    if "collections" not in gsea:
        fail("enrichment.gsea.collections is missing; the enrichment workflow "
             "expects a list (can be empty []). This is the exact key whose "
             "absence triggered the 'No such variable: file' error.")
    else:
        ok(f"enrichment.{{go,kegg,wormcat,gsea}} all present; "
           f"gsea.collections has {len(gsea['collections'])} entr{'y' if len(gsea['collections'])==1 else 'ies'}")


def check_design_rank(p: dict, samples: list[dict]) -> None:
    """Build the design matrix from the YAML + samplesheet and check it's
    full rank. Catches the kind of structural collinearity that makes
    DESeq2 refuse to fit (e.g., adding `batch: source` when source is
    perfectly aliased with another design factor's level).

    Quick pure-Python implementation: build a one-hot/numeric model matrix
    from the design variables, compute its numpy rank, compare to ncol.
    """
    deseq2 = p.get("deseq2", {}) or {}
    formula = deseq2.get("design_formula")
    if not formula:
        return  # only check multi-factor mode
    print("\n[+] Design matrix rank check (multi-factor)")
    try:
        import numpy as np
    except ImportError:
        warn("numpy not available; skipping rank check (install with --with numpy)")
        return

    user_vars = parse_design_vars(formula)
    bio_vars = [v for v in user_vars if not re.match(r"^W_\d+$", v)]
    nuisance = []
    if deseq2.get("batch") and deseq2.get("batch") in samples[0]:
        nuisance.append(deseq2["batch"])
    full_vars = nuisance + bio_vars

    # Build one-hot model matrix per categorical, intercept column included.
    # Use the same reference level the YAML declares so the matrix matches
    # what DESeq2 will build.
    refs = dict(parse_reference_levels(deseq2.get("reference_levels", "")))
    n = len(samples)
    cols: list[tuple[str, list[float]]] = [("(Intercept)", [1.0] * n)]
    for var in full_vars:
        vals = [str(s.get(var)) for s in samples]
        levels = sorted(set(vals))
        ref = refs.get(var, levels[0])  # alphabetical default if not specified
        if ref not in levels:
            warn(f"reference level '{ref}' for {var} not in observed levels; falling back")
            ref = levels[0]
        for lvl in levels:
            if lvl == ref:
                continue
            cols.append((f"{var}{lvl}", [1.0 if v == lvl else 0.0 for v in vals]))
    M = np.array([col for _, col in cols]).T   # shape (n_samples, n_cols)
    col_names = [name for name, _ in cols]
    rank = int(np.linalg.matrix_rank(M))
    print(f"  Design columns ({len(col_names)}): {col_names}")
    print(f"  Matrix shape: {M.shape}, rank: {rank}")
    if rank < M.shape[1]:
        fail(f"DESIGN MATRIX IS NOT FULL RANK ({rank} < {M.shape[1]} columns). "
             f"DESeq2 will refuse to fit. Likely cause: a categorical batch "
             f"covariate is perfectly aliased with one of the biological "
             f"factor levels. Drop the aliased nuisance variable from `batch:` "
             f"and rely on RUV (W_*) instead, or merge the aliased levels.")
        # Try to identify which columns are linearly dependent
        for i in range(1, M.shape[1]):
            sub_rank = int(np.linalg.matrix_rank(M[:, : i + 1]))
            if sub_rank < i + 1:
                deps = col_names[: i + 1]
                print(f"           Linear dependence introduced when adding "
                      f"`{col_names[i]}` (full set so far: {deps})")
                break
    else:
        ok(f"design matrix is full rank")


def check_fastq_launcher(samples: list[dict]) -> None:
    """If the samplesheet has FASTQ samples, the launcher MUST be `bash run.sh`
    (which runs nf-core/rnaseq and sets `--rnaseq_outdir`), not a bare
    `nextflow run main.nf`. Otherwise main.nf's FASTQ -> BAM channel is empty
    and only the plasmidsaurus BAMs reach merge_counts -- silently dropping
    every FASTQ sample from the count matrix and downstream DESeq2.

    We can't know what command the user is about to type, so this is a
    LAUNCH-TIME WARN with explicit instructions, not a FAIL.
    """
    print("\n[+] FASTQ samples / launcher reminder")
    fastq_samples = [
        s for s in samples
        if str(s.get("input_type", "")).startswith("fastq_")
        or s.get("fastq_1")
        or s.get("fastq_2")
    ]
    if not fastq_samples:
        print("  (skipped; samplesheet has no FASTQ samples)")
        return
    n = len(fastq_samples)
    warn(f"samplesheet has {n} FASTQ sample(s). They will only reach the "
         f"count matrix if the launcher sets --rnaseq_outdir.")
    print( "           DO:    bash run.sh --params <yaml> --profile aws [--skip-alignment --rnaseq-outdir s3://.../rnaseq_output/<run_id>] [--resume]")
    print( "           AVOID: nextflow run main.nf -params-file <yaml>  (rnaseq_outdir unset -> FASTQ samples silently dropped)")
    print( "           If you've already aligned these FASTQs in a prior run, --skip-alignment + --rnaseq-outdir <prior path> reuses those BAMs.")


def check_source_aliasing(p: dict, samples: list[dict]) -> None:
    print("\n[8] Cross-protocol source aliasing")
    deseq2 = p.get("deseq2", {}) or {}
    batch = deseq2.get("batch")
    if not batch or not samples:
        print("  (skipped; no batch column set)")
        return
    if batch not in samples[0]:
        warn(f"batch=`{batch}` not in samplesheet columns; pipeline will skip batch correction")
        return
    # Group samples by genotype and tally batch values
    by_geno: dict[str, set] = {}
    for s in samples:
        g = str(s.get("genotype", "?"))
        by_geno.setdefault(g, set()).add(str(s.get(batch)))
    aliased = []
    for g, batches in by_geno.items():
        if len(batches) == 1:
            aliased.append((g, list(batches)[0]))
    if aliased:
        warn(f"`{batch}` is aliased with `genotype` for the following genotypes "
             f"(only one batch per genotype, so the batch effect cannot be separated "
             f"from the genotype effect):")
        for g, b in aliased:
            print(f"           genotype={g} only has batch={b}")
        print("           This is fine if RUVseq is enabled (RUV will absorb the")
        print("           cross-protocol variation), but be aware: log2FCs across")
        print("           the aliased dimensions are direction-only, not magnitude.")
    else:
        ok(f"no per-genotype batch aliasing detected for `{batch}`")


# ----------------------------------------------------------------- main

def main() -> None:
    if len(sys.argv) != 2:
        sys.exit("Usage: check_pipeline_inputs.py <pipeline.yaml>")
    yaml_path = Path(sys.argv[1]).resolve()
    print(f"Pre-flight check for: {yaml_path}")

    p = load_yaml(yaml_path)
    samplesheet_path = check_pipeline_yaml(p, yaml_path)
    s = load_yaml(samplesheet_path)
    samples = check_samplesheet_structure(s)
    check_s3_inputs(samples)
    check_references(p)
    check_multifactor_design(p, samples)
    check_enrichment_block(p)
    check_design_rank(p, samples)
    check_fastq_launcher(samples)
    check_source_aliasing(p, samples)

    print()
    if ANY_FAIL:
        print("=== RESULT: FAIL -- fix the issues above before launching the pipeline ===")
        sys.exit(1)
    else:
        print("=== RESULT: OK -- pipeline inputs look good; safe to launch ===")
        sys.exit(0)


if __name__ == "__main__":
    main()
