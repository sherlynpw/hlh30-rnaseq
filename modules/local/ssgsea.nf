process SSGSEA {
    tag "ssgsea"
    label 'process_medium'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/enrichment/ssgsea", mode: params.publish_dir_mode

    input:
    path vst_counts
    path sample_meta
    path gene_db
    path manifest_json
    path collection_files
    path ortholog_map
    path de_csvs, stageAs: 'de_csvs/*'   // flat dir of DESeq2 CSVs (default + compound contrasts)

    output:
    path "**/*.csv", emit: tables, optional: true
    path "**/*.pdf", emit: plots,  optional: true
    path "**",       emit: results, optional: true

    when:
    params.enrichment.ssgsea.enabled

    script:
    def cfg = params.enrichment.ssgsea
    def gsea_cfg = params.enrichment.gsea
    // Reuse GSEA builtin list (go/kegg/...) unless ssgsea.builtin overrides.
    // gsea.collections is a list of external-GMT maps and is handled via the manifest.
    def builtin_src = cfg.builtin != null ? cfg.builtin : (gsea_cfg.builtin ?: [])
    def builtin_list = (builtin_src ?: []).join(',')
    def builtin_arg = builtin_list ? "--builtin ${builtin_list}" : ''
    def manifest_arg = manifest_json.name != 'NO_COLLECTIONS_MANIFEST' ? "--collections_manifest ${manifest_json}" : ''
    def db_arg = gene_db.name != 'NO_GENE_DB' ? "--gene_db ${gene_db}" : ''
    def ortholog_map_arg = ortholog_map.name != 'NO_ORTHOLOG_MAP' ? "--ortholog_map ${ortholog_map}" : ''
    def ortholog_from_arg = gsea_cfg.ortholog_from ? "--ortholog_from ${gsea_cfg.ortholog_from}" : ''
    def ortholog_to_arg   = gsea_cfg.ortholog_to   ? "--ortholog_to ${gsea_cfg.ortholog_to}"     : ''
    def min_size_arg = cfg.min_size ? "--min_size ${cfg.min_size}" : ''
    def max_size_arg = cfg.max_size ? "--max_size ${cfg.max_size}" : ''
    def top_n_arg    = cfg.heatmap_top_n ? "--top_n ${cfg.heatmap_top_n}" : ''
    def method_arg   = cfg.method ? "--method ${cfg.method}" : ''
    def formula_arg  = cfg.formula ? "--formula '${cfg.formula}'" : ''
    def coef_list    = (cfg.coefficients ?: []).join(',')
    def coef_arg     = coef_list ? "--coefficients '${coef_list}'" : ''
    def color_list   = (cfg.color_by ?: []).join(',')
    def color_arg    = color_list ? "--color_by '${color_list}'" : ''
    def ref_list     = (cfg.reference_levels ?: [:]).collect { k, v -> "${k}=${v}" }.join(',')
    def ref_arg      = ref_list ? "--reference_levels '${ref_list}'" : ''
    def top_paths_arg = cfg.top_pathways ? "--top_pathways ${cfg.top_pathways}" : ''
    def dp_top_arg    = cfg.dotplot_top_n    ? "--dotplot_top_n ${cfg.dotplot_top_n}"       : ''
    def vol_lab_arg   = cfg.volcano_label_n  ? "--volcano_label_n ${cfg.volcano_label_n}"   : ''
    def traj_x_arg   = cfg.trajectory_x ? "--trajectory_x ${cfg.trajectory_x}" : ''
    def traj_grp_arg = cfg.trajectory_group ? "--trajectory_group ${cfg.trajectory_group}" : ''
    def box_x_list   = (cfg.box_x ?: []).join(',')
    def box_x_arg    = box_x_list ? "--box_x '${box_x_list}'" : ''
    // Column ordering: default to clustering.condition_order over (genotype, timepoint)-style cols
    def order_cols_list = (cfg.column_order_cols ?: []).join(',')
    def order_cols_arg  = order_cols_list ? "--column_order_cols '${order_cols_list}'" : ''
    def order_csv = cfg.column_order != null ? (cfg.column_order instanceof List ? cfg.column_order.join(',') : cfg.column_order)
                                              : (params.clustering?.condition_order ?: '')
    def order_arg = order_csv ? "--column_order '${order_csv}'" : ''
    // Named linear contrasts → JSON file the R script reads
    def contrasts_map = cfg.contrasts ?: [:]
    def contrasts_json = contrasts_map ? groovy.json.JsonOutput.toJson(contrasts_map) : ''
    def contrasts_setup = contrasts_json ? "cat > contrasts.json <<'_EOF_'\n${contrasts_json}\n_EOF_" : ''
    def contrasts_arg = contrasts_json ? "--contrasts_json contrasts.json" : ''
    // DE mapping: ssGSEA term name -> DE CSV basename
    def de_mapping = cfg.de_mapping ?: [:]
    def de_mapping_json = de_mapping ? groovy.json.JsonOutput.toJson(de_mapping) : ''
    def de_mapping_setup = de_mapping_json ? "cat > de_mapping.json <<'_EOF_'\n${de_mapping_json}\n_EOF_" : ''
    def de_mapping_arg = de_mapping_json ? "--de_mapping_json de_mapping.json" : ''
    def de_dir_arg = de_csvs ? "--de_results_dir de_csvs" : ''
    def driver_top_arg = cfg.driver_top_n ? "--driver_top_n ${cfg.driver_top_n}" : ''
    """
    ${contrasts_setup}
    ${de_mapping_setup}
    ssgsea.R \\
        --vst_counts ${vst_counts} \\
        --sample_meta ${sample_meta} \\
        --org_db ${params.organism.org_db} \\
        --kegg_org ${params.organism.kegg_org ?: 'cel'} \\
        ${builtin_arg} \\
        ${manifest_arg} \\
        ${db_arg} \\
        ${ortholog_map_arg} \\
        ${ortholog_from_arg} \\
        ${ortholog_to_arg} \\
        ${min_size_arg} \\
        ${max_size_arg} \\
        ${top_n_arg} \\
        ${method_arg} \\
        ${formula_arg} \\
        ${coef_arg} \\
        ${color_arg} \\
        ${ref_arg} \\
        ${top_paths_arg} \\
        ${dp_top_arg} \\
        ${vol_lab_arg} \\
        ${traj_x_arg} \\
        ${traj_grp_arg} \\
        ${box_x_arg} \\
        ${order_cols_arg} \\
        ${order_arg} \\
        ${contrasts_arg} \\
        ${de_dir_arg} \\
        ${de_mapping_arg} \\
        ${driver_top_arg} \\
        --outdir .
    """
}
