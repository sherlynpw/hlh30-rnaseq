FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL maintainer="Ray Andrew"
LABEL description="C. elegans RNA-seq pipeline R analysis container"

# Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'DESeq2', \
    'clusterProfiler', \
    'org.Ce.eg.db', \
    'biomaRt', \
    'Mfuzz', \
    'sva', \
    'limma', \
    'RUVSeq', \
    'EDASeq', \
    'fgsea', \
    'GSVA', \
    'enrichplot', \
    'AnnotationDbi', \
    'GO.db', \
    'Rsubread', \
    'apeglm', \
    'ashr', \
    'motifmatchr', \
    'TFBSTools', \
    'JASPAR2024', \
    'BSgenome.Celegans.UCSC.ce11', \
    'GenomicFeatures', \
    'txdbmaker', \
    'rtracklayer', \
    'Biostrings', \
    'Gviz', \
    'ComplexHeatmap', \
    'monaLisa', \
    'rrvgo' \
), ask=FALSE, update=FALSE, force=TRUE)"

# CRAN packages
RUN R -e "install.packages(c( \
    'optparse', \
    'dplyr', \
    'readr', \
    'tidyr', \
    'ggplot2', \
    'pheatmap', \
    'patchwork', \
    'uwot', \
    'dbscan', \
    'dynamicTreeCut', \
    'RSQLite', \
    'jsonlite', \
    'yaml', \
    'scales', \
    'tibble', \
    'ggseqlogo', \
    'ggvenn', \
    'ggrepel', \
    'treemap', \
    'wordcloud' \
), repos='https://cloud.r-project.org')"

# wormcat from GitHub + pre-download annotation files
RUN R -e "library(devtools); install_github('dphiggs01/wormcat')" && \
    R -e "library(wormcat); copy_annotation_files_to_extdata()"

# Verify ALL packages used by pipeline scripts
RUN R -e " \
    pkgs <- c( \
        'DESeq2', 'clusterProfiler', 'org.Ce.eg.db', 'biomaRt', \
        'Mfuzz', 'Biobase', 'fgsea', 'GSVA', 'uwot', 'dbscan', 'sva', 'limma', \
        'RUVSeq', 'EDASeq', \
        'enrichplot', 'AnnotationDbi', 'GO.db', 'apeglm', 'Rsubread', \
        'optparse', 'dplyr', 'tidyr', 'readr', 'ggplot2', 'pheatmap', 'patchwork', \
        'dynamicTreeCut', \
        'RSQLite', 'jsonlite', 'yaml', 'scales', 'tibble', 'wormcat', 'ashr', \
        'motifmatchr', 'TFBSTools', 'JASPAR2024', 'BSgenome.Celegans.UCSC.ce11', \
        'GenomicFeatures', 'txdbmaker', 'rtracklayer', 'Biostrings', 'Gviz', \
        'ComplexHeatmap', 'monaLisa', 'ggseqlogo', 'ggvenn', 'ggrepel', \
        'rrvgo', 'treemap', 'wordcloud' \
    ); \
    for (p in pkgs) { \
        if (!requireNamespace(p, quietly=TRUE)) stop(paste('MISSING:', p)); \
    }; \
    message('All ', length(pkgs), ' packages verified') \
"

# AWS CLI (for S3 staging on AWS Batch)
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip -q awscliv2.zip && ./aws/install && \
    rm -rf awscliv2.zip aws/

WORKDIR /data
