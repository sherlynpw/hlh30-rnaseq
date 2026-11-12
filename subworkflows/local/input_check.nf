import groovy.yaml.YamlSlurper

workflow INPUT_CHECK {

    take:
    samplesheet  // path to YAML samplesheet

    main:
    // Parse YAML samplesheet
    Channel
        .fromPath(samplesheet)
        .map { yaml_file ->
            def yaml = new YamlSlurper().parse(yaml_file)
            yaml.samples
        }
        .flatMap { it }
        .map { sample ->
            // Build meta from all samplesheet columns (except file paths)
            def skip_keys = ['bam', 'fastq_1', 'fastq_2'] as Set
            def meta = [id: sample.sample_id, input_type: sample.input_type]
            sample.each { k, v ->
                if (k != 'sample_id' && k != 'input_type' && !(k in skip_keys)) {
                    meta[k] = v
                }
            }
            // Defaults
            meta.source       = meta.source       ?: 'unknown'
            meta.strandedness = meta.strandedness  ?: 'auto'
            [meta, sample]
        }
        .branch { meta, sample ->
            plasmidsaurus_bam: meta.input_type == 'plasmidsaurus_bam'
                return [meta, file(sample.bam)]
            fastq_3prime: meta.input_type == 'fastq_3prime'
                return [meta, sample.fastq_2 ? [file(sample.fastq_1), file(sample.fastq_2)] : [file(sample.fastq_1)]]
            fastq_standard: meta.input_type == 'fastq_standard'
                return [meta, sample.fastq_2 ? [file(sample.fastq_1), file(sample.fastq_2)] : [file(sample.fastq_1)]]
        }
        .set { ch_branched }

    emit:
    plasmidsaurus_bams = ch_branched.plasmidsaurus_bam  // [meta, bam]
    fastq_3prime       = ch_branched.fastq_3prime        // [meta, [fastq]]
    fastq_standard     = ch_branched.fastq_standard      // [meta, [fastq]]
}
