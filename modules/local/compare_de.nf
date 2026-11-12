process COMPARE_DE {
    tag "compare_de"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/comparison", mode: params.publish_dir_mode

    input:
    path de_results

    output:
    path "compare_de_report.csv",  emit: report
    path "*.pdf",                  emit: plots, optional: true
    path "*",                      emit: results

    script:
    // Optional age/genotype ordering. If omitted, compare_de.R auto-derives
    // alphabetically. For workloads where alphabetical doesn't match biology
    // (e.g. L4 < Day3 < Day6), set enrichment.compare_de_age_order in the
    // params YAML, comma-separated, oldest LAST.
    def age_arg   = params.enrichment?.compare_de_age_order  ? "--age_order '${params.enrichment.compare_de_age_order}'"   : ''
    def geno_arg  = params.enrichment?.compare_de_geno_order ? "--geno_order '${params.enrichment.compare_de_geno_order}'" : ''

    // Optional pair allow-list: enrichment.compare_de_pairs is a list of
    // "A_vs_B__X__C_vs_D" strings. Only those cross-comparisons emit results.
    def pairs_arg = params.enrichment?.compare_de_pairs ?
        "--pairs '${params.enrichment.compare_de_pairs.join(',')}'" : ''

    // Optional term allow-list: enrichment.compare_de_include_terms is a list
    // of "<genotype>_<age>" tokens. Result files whose num and den are both
    // in the list are kept; others are skipped before the pair iteration.
    // In canonicalize=false mode, this becomes a factor-name allow-list
    // (e.g. ["genotype"] to drop timepoint/replicate coefs).
    def terms_arg = params.enrichment?.compare_de_include_terms ?
        "--include_terms '${params.enrichment.compare_de_include_terms.join(',')}'" : ''

    // Canonicalize flag. Default TRUE (legacy genotype/age parsing). Set to
    // FALSE in multi-factor pipeline configs where contrast files are named
    // <factor>__<level>_vs_<reference> and don't have a genotype/age axis.
    def canon_arg = (params.enrichment?.compare_de_canonicalize == false) ?
        "--canonicalize FALSE" : "--canonicalize TRUE"

    // Comparison method: 'overlap' (legacy classification) or 'ztest' (Wald
    // test on the LFC difference). Defaults to overlap when unset.
    def method_arg = params.enrichment?.compare_de_method ?
        "--diff_method ${params.enrichment.compare_de_method}" : ''
    """
    compare_de.R \\
        --input_dir . \\
        --outdir . \\
        ${age_arg} \\
        ${geno_arg} \\
        ${terms_arg} \\
        ${pairs_arg} \\
        ${canon_arg} \\
        ${method_arg}
    """
}
