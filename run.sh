#!/usr/bin/env bash
set -euo pipefail

# C. elegans RNA-seq Pipeline - Wrapper Script
#
# Orchestrates nf-core/rnaseq (for FASTQ inputs) + custom analysis pipeline.
# For Plasmidsaurus BAM-only runs, skip step 1 entirely.
#
# Default behavior:
#   --profile aws   -> --rnaseq-outdir defaults to s3://${AWS_BUCKET}/rnaseq_output
#                      (BAMs/QC stay on S3; main.nf reads them straight from S3)
#   --profile local -> --rnaseq-outdir defaults to ./rnaseq_output (local dir)
#
# A small local --staging-dir (default ./rnaseq_staging) holds the generated
# nf-core CSV samplesheet + Nextflow's local CWD/.nextflow history; both must
# be local even when outdir is on S3.
#
# Usage:
#   ./run.sh --params params.yaml --profile local
#   ./run.sh --params params.yaml --profile aws
#   ./run.sh --params params.yaml --profile aws --skip-alignment   # reuse existing S3 BAMs
#   ./run.sh --params params.yaml --profile aws --rnaseq-outdir s3://other/path
#   ./run.sh --params params.yaml --profile aws --resume

PARAMS_FILE="params.yaml"
PROFILE="local"
SKIP_ALIGNMENT=false
RNASEQ_OUTDIR=""           # auto-defaults below based on profile
LOCAL_STAGING_DIR="./rnaseq_staging"
RESUME=""
DRY_RUN=false
NFCORE_VERSION="3.14.0"

while [[ $# -gt 0 ]]; do
    case $1 in
        --params)         PARAMS_FILE="$2"; shift 2 ;;
        --profile)        PROFILE="$2"; shift 2 ;;
        --skip-alignment) SKIP_ALIGNMENT=true; shift ;;
        --rnaseq-outdir)  RNASEQ_OUTDIR="$2"; shift 2 ;;
        --staging-dir)    LOCAL_STAGING_DIR="$2"; shift 2 ;;
        --resume)         RESUME="-resume"; shift ;;
        --dry-run)        DRY_RUN=true; shift ;;
        *)                echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Reference paths - read from env.yaml so nf-core/rnaseq aligns to the SAME
# Ensembl release that main.nf counts against (otherwise BAM contigs won't
# match the GTF and featureCounts silently returns ~0 reads).
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/scripts/env.sh"

REF_GTF="${REF_GTF:-s3://${AWS_BUCKET}/refs/Caenorhabditis_elegans.${GENOME}.${ENSEMBL_RELEASE}.gtf}"
REF_FASTA="${REF_FASTA:-s3://${AWS_BUCKET}/refs/Caenorhabditis_elegans.${GENOME}.dna.toplevel.fa}"
STAR_INDEX="${STAR_INDEX:-s3://${AWS_BUCKET}/refs/star_${GENOME}_${ENSEMBL_RELEASE}/}"   # set STAR_INDEX= to disable

# Pull `id:` out of the params YAML. Used only for local staging dirs so
# concurrent runs of different pipelines don't clobber each other's generated
# CSVs / .nextflow history. BAM output is shared across runs (pool layout),
# protected from cross-run sample leakage by main.nf's samplesheet filter.
RUN_ID="$(grep -E '^id:' "${PARAMS_FILE}" | head -1 | awk '{print $2}' | tr -d '"')"
if [[ -z "${RUN_ID}" ]]; then
    echo "ERROR: could not parse 'id:' from ${PARAMS_FILE}" >&2
    exit 1
fi

# Default rnaseq-outdir: shared BAM pool. Same SRR aligned across multiple
# pipelines lives in one place so the second pipeline doesn't re-align it.
# Pass --rnaseq-outdir to override (e.g. for an isolated experiment).
if [[ -z "${RNASEQ_OUTDIR}" ]]; then
    case "${PROFILE}" in
        aws) RNASEQ_OUTDIR="s3://${AWS_BUCKET}/rnaseq_output/_pool" ;;
        *)   RNASEQ_OUTDIR="./rnaseq_output/_pool" ;;
    esac
fi

# Local staging is also per-id so concurrent runs don't clobber each other's
# generated samplesheets / .nfcore-cwd
LOCAL_STAGING_DIR="${LOCAL_STAGING_DIR%/}/${RUN_ID}"

echo "=== C. elegans RNA-seq Pipeline ==="
echo "Params:         ${PARAMS_FILE}"
echo "Run id:         ${RUN_ID}"
echo "Profile:        ${PROFILE}"
echo "Skip alignment: ${SKIP_ALIGNMENT}"
echo "Ref FASTA:      ${REF_FASTA}"
echo "Ref GTF:        ${REF_GTF}"
[[ -n "${STAR_INDEX}" ]] && echo "STAR index:     ${STAR_INDEX}"
echo "Rnaseq outdir:  ${RNASEQ_OUTDIR}"
echo "Local staging:  ${LOCAL_STAGING_DIR}"
echo ""

# Step 1: Run nf-core/rnaseq for any FASTQ samples
if [[ "${SKIP_ALIGNMENT}" == false ]]; then
    # Check if there are FASTQ samples in the samplesheet
    # This is a simple grep - the pipeline itself handles proper routing
    if grep -q 'fastq_' "${PARAMS_FILE}" 2>/dev/null || \
       grep -q 'fastq_' "$(grep 'input:' "${PARAMS_FILE}" | awk '{print $2}')" 2>/dev/null; then
        echo "=== Step 1: Running nf-core/rnaseq for FASTQ samples ==="

        # Samplesheet CSVs and the Nextflow CWD MUST be local (Nextflow needs
        # writable local CWD; --input expects a local CSV). RNASEQ_OUTDIR can
        # be S3 - Nextflow publishes there, and main.nf reads BAMs from it.
        mkdir -p "${LOCAL_STAGING_DIR}"
        Rscript bin/generate_nfcore_samplesheet.R \
            --input "$(grep 'input:' "${PARAMS_FILE}" | awk '{print $2}')" \
            --outdir "${LOCAL_STAGING_DIR}"

        # Auto-skip alignment if every expected BAM is already published.
        # Reads sample names from the freshly-generated nf-core CSV(s) and
        # checks the corresponding *.markdup.sorted.bam at <outdir>/{standard,
        # 3prime}/star_salmon/. If all present, skip nf-core/rnaseq entirely.
        check_bams_present() {
            local sheet="$1" subdir="$2"
            [[ -f "$sheet" ]] || return 0
            local samples
            samples=$(awk -F, 'NR>1 && $1!=""{print $1}' "$sheet")
            [[ -z "$samples" ]] && return 0
            local missing=0
            local bam_root="${RNASEQ_OUTDIR}/${subdir}/star_salmon"
            for s in $samples; do
                local bam="${bam_root}/${s}.markdup.sorted.bam"
                if [[ "$bam" == s3://* ]]; then
                    aws s3 ls "$bam" >/dev/null 2>&1 || { missing=1; break; }
                else
                    [[ -f "$bam" ]] || { missing=1; break; }
                fi
            done
            return $missing
        }
        STD_SHEET_PROBE="${LOCAL_STAGING_DIR}/samplesheet_standard.csv"
        TPRIME_SHEET_PROBE="${LOCAL_STAGING_DIR}/samplesheet_3prime.csv"
        all_present=true
        if [[ -f "${STD_SHEET_PROBE}" ]] && ! check_bams_present "${STD_SHEET_PROBE}" "standard"; then
            all_present=false
        fi
        if [[ -f "${TPRIME_SHEET_PROBE}" ]] && ! check_bams_present "${TPRIME_SHEET_PROBE}" "3prime"; then
            all_present=false
        fi
        if [[ "$all_present" == true ]]; then
            echo "All expected BAMs already exist under ${RNASEQ_OUTDIR} - auto-skipping alignment."
            SKIP_ALIGNMENT=true
        else
            # Pool layout: drop CSV rows whose BAM already lives in the pool so
            # nf-core only re-aligns the missing samples. The full row set is
            # still validated by check_bams_present after nf-core completes via
            # main.nf's samplesheet filter, so a dropped row never silently
            # disappears from the analysis.
            filter_csv_to_missing() {
                local sheet="$1" subdir="$2"
                [[ -f "$sheet" ]] || return 0
                local bam_root="${RNASEQ_OUTDIR}/${subdir}/star_salmon"
                local tmp="${sheet}.filtered"
                head -1 "$sheet" > "$tmp"
                local kept=0 skipped=0
                while IFS= read -r line; do
                    local s
                    s=$(echo "$line" | awk -F, '{print $1}')
                    [[ -z "$s" ]] && continue
                    local bam="${bam_root}/${s}.markdup.sorted.bam"
                    local present=1
                    if [[ "$bam" == s3://* ]]; then
                        aws s3 ls "$bam" >/dev/null 2>&1 || present=0
                    else
                        [[ -f "$bam" ]] || present=0
                    fi
                    if [[ "$present" == 1 ]]; then
                        skipped=$((skipped+1))
                    else
                        echo "$line" >> "$tmp"
                        kept=$((kept+1))
                    fi
                done < <(tail -n +2 "$sheet")
                if [[ "$skipped" -gt 0 ]]; then
                    echo "  ${subdir}: ${skipped} BAM(s) already in pool, ${kept} to align"
                    mv "$tmp" "$sheet"
                else
                    rm -f "$tmp"
                fi
            }
            echo "Checking pool for partial reuse:"
            filter_csv_to_missing "${STD_SHEET_PROBE}" "standard"
            filter_csv_to_missing "${TPRIME_SHEET_PROBE}" "3prime"
        fi
    fi
fi

if [[ "${DRY_RUN}" == true ]]; then
    echo ""
    echo "=== Dry run ==="
    echo "Skip alignment:  ${SKIP_ALIGNMENT}"
    echo "Pool root:       ${RNASEQ_OUTDIR}"
    for sheet in "${LOCAL_STAGING_DIR}/samplesheet_standard.csv" \
                 "${LOCAL_STAGING_DIR}/samplesheet_3prime.csv"; do
        if [[ -f "$sheet" ]]; then
            nrows=$(( $(wc -l < "$sheet") - 1 ))
            echo ""
            echo "  ${sheet}:  ${nrows} sample(s) would be aligned by nf-core"
            [[ "$nrows" -gt 0 ]] && tail -n +2 "$sheet" | awk -F, '{print "    - "$1}'
        fi
    done
    echo ""
    echo "main.nf would then read pool BAMs from ${RNASEQ_OUTDIR}/{standard,3prime}/star_salmon/"
    echo "(filtered to sample_ids declared in $(grep '^input:' "${PARAMS_FILE}" | awk '{print $2}')) "
    echo ""
    echo "No alignment or analysis was executed. Drop --dry-run to actually run."
    exit 0
fi

# Actually run nf-core/rnaseq if alignment is still needed.
if [[ "${SKIP_ALIGNMENT}" == false ]]; then
    STAR_INDEX_ARG=""
    [[ -n "${STAR_INDEX}" ]] && STAR_INDEX_ARG="--star_index ${STAR_INDEX}"

    # Avoid leaking our project's nextflow.config into nf-core/rnaseq.
    # Nextflow auto-loads nextflow.config from CWD; ours declares our
    # pipeline params (id, organism, ...). nf-validation in nf-core/rnaseq
    # rejects unknown params -> "Cannot get property 'properties' on null".
    # Solution: run nf-core from a clean local CWD with our infra-only conf.
    NFCORE_DIR="${LOCAL_STAGING_DIR}/.nfcore-cwd"
    mkdir -p "${NFCORE_DIR}"
    case "${PROFILE}" in
        aws)   NFCORE_EXTRA_CONF=( -c "${SCRIPT_DIR}/conf/nfcore-rnaseq-aws.config" ) ;;
        local) NFCORE_EXTRA_CONF=( -profile docker ) ;;
        *)     NFCORE_EXTRA_CONF=( -profile "${PROFILE}" ) ;;
    esac

    STD_SHEET="$(realpath "${LOCAL_STAGING_DIR}/samplesheet_standard.csv" 2>/dev/null || echo "")"
    TPRIME_SHEET="$(realpath "${LOCAL_STAGING_DIR}/samplesheet_3prime.csv" 2>/dev/null || echo "")"

    if [[ "${RNASEQ_OUTDIR}" == s3://* ]]; then
        STD_OUT="${RNASEQ_OUTDIR}/standard"
        TPRIME_OUT="${RNASEQ_OUTDIR}/3prime"
    else
        STD_OUT="$(realpath -m "${RNASEQ_OUTDIR}/standard")"
        TPRIME_OUT="$(realpath -m "${RNASEQ_OUTDIR}/3prime")"
    fi

    # Header-only CSV (all samples already in pool for this input type) -> skip
    csv_has_rows() { [[ -f "$1" ]] && [[ "$(wc -l < "$1")" -ge 2 ]]; }

    if [[ -n "${STD_SHEET}" ]] && csv_has_rows "${STD_SHEET}"; then
        ( cd "${NFCORE_DIR}" && \
          nextflow run nf-core/rnaseq \
            -r "${NFCORE_VERSION}" \
            "${NFCORE_EXTRA_CONF[@]}" \
            --input "${STD_SHEET}" \
            --outdir "${STD_OUT}" \
            --aligner star_salmon \
            --fasta "${REF_FASTA}" \
            --gtf "${REF_GTF}" \
            ${STAR_INDEX_ARG} \
            --save_align_intermeds \
            --skip_qc \
            --skip_multiqc \
            --skip_pseudo_alignment \
            ${RESUME} )
    fi

    if [[ -n "${TPRIME_SHEET}" ]] && csv_has_rows "${TPRIME_SHEET}"; then
        ( cd "${NFCORE_DIR}" && \
          nextflow run nf-core/rnaseq \
            -r "${NFCORE_VERSION}" \
            "${NFCORE_EXTRA_CONF[@]}" \
            --input "${TPRIME_SHEET}" \
            --outdir "${TPRIME_OUT}" \
            --aligner star_salmon \
            --fasta "${REF_FASTA}" \
            --gtf "${REF_GTF}" \
            ${STAR_INDEX_ARG} \
            --save_align_intermeds \
            --skip_qc \
            --skip_multiqc \
            --skip_pseudo_alignment \
            --trimmer fastp \
            --extra_fastp_args '--trim_poly_x --poly_x_min_len 10' \
            ${RESUME} )
    fi

    echo "=== Step 1 complete ==="
fi

# Step 2: Run custom analysis pipeline.
# Run from a per-id CWD so each pipeline gets its own .nextflow/ history and
# session UUID -- otherwise two parallel `--resume` runs of main.nf from the
# project root collide on the same .nextflow/cache lock.
echo "=== Step 2: Running custom analysis pipeline ==="
MAIN_DIR="${LOCAL_STAGING_DIR}/.main-cwd"
mkdir -p "${MAIN_DIR}"
MAIN_NF="$(realpath "${SCRIPT_DIR}/main.nf")"
PARAMS_ABS="$(realpath "${PARAMS_FILE}")"

# Pre-resolve `input:` path from the params YAML to absolute (since main.nf
# resolves it relative to CWD, and we just changed CWD to .main-cwd).
INPUT_REL="$(grep '^input:' "${PARAMS_FILE}" | head -1 | awk '{print $2}' | tr -d '"')"
INPUT_ABS="$(realpath "${SCRIPT_DIR}/${INPUT_REL}")"

( cd "${MAIN_DIR}" && \
  nextflow run "${MAIN_NF}" \
    -params-file "${PARAMS_ABS}" \
    --input "${INPUT_ABS}" \
    --rnaseq_outdir "${RNASEQ_OUTDIR}" \
    -profile "${PROFILE}" \
    ${RESUME} )

echo "=== Pipeline complete ==="
