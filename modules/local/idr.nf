process idr {
    tag "$meta.id"
    label 'process_high'

    // custom IDR measurement/filtering of peaks grouped by antibody (as listed in input csv). Adapted from bioinfo-pf-curie/ChIP-seq and ENCODE ChIP-seq, which use https://github.com/nboley/idr

    conda (params.enable_conda ? "bioconda::idr=2.0.4.2=py39h5371cbf_7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/idr:2.0.4.2--py39h5371cbf_7' :
        'quay.io/biocontainers/idr:2.0.4.2--py39h5371cbf_7' }"

    input:
    tuple val(meta), path(combined_peaks)
    val idr_cutoff
    
    output:
    path "*.bed"         , emit: idr 
    path "*log.txt"      , emit: mqcIdr 
    path "v_idr.txt"     , emit: version
    path "*.png"         , emit: plot
    
    script:
    peaktype = combined_peaks[0].toString()
    peaktype = peaktype.substring(peaktype.lastIndexOf(".") + 1)
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ 
    idr --version &> v_idr.txt
    idr --samples ${combined_peaks} \\
        --input-file-type bed \\
        --rank 5 \\
        --idr-threshold ${idr_cutoff} \\
        --output-file-type bed \\
        -o ${prefix}.bed \\
        -l ${prefix}_log.txt \\
        --plot
    """
}
