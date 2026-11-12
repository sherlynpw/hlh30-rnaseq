include { PCA                } from '../../modules/local/pca'
include { UMAP               } from '../../modules/local/umap'
include { CV_DENSITY         } from '../../modules/local/cv_density'
include { SAMPLE_CORRELATION } from '../../modules/local/sample_correlation'

workflow EXPLORATION {

    take:
    ch_vst       // VST-normalized counts (corrected if clustering.correct_vst)
    ch_vst_raw   // Always-raw VST for before/after diagnostics
    ch_meta      // sample metadata

    main:
    PCA(ch_vst, ch_vst_raw, ch_meta)
    UMAP(ch_vst, ch_vst_raw, ch_meta)
    CV_DENSITY(ch_vst)
    SAMPLE_CORRELATION(ch_vst, ch_meta)

    emit:
    pca_results        = PCA.out.results
    umap_results       = UMAP.out.results
    cv_density         = CV_DENSITY.out.plots
    sample_correlation = SAMPLE_CORRELATION.out.results
}
