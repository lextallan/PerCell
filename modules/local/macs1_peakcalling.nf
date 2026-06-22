process macs1_peakcalling {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::macs==1.4.3--pyhdfd78af_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macs:1.4.3--pyhdfd78af_0' :
        'quay.io/biocontainers/macs:1.4.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ipbam), path(controlbam)
    val macs_gsize
    val macs1_pvalue

    output:
    tuple val(meta), path("*peaks.bed"), emit: peak
    tuple val(meta), path("*summits.bed"), emit: summit

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    macs -t ${ipbam} -c ${controlbam} -p ${macs1_pvalue} -g ${macs_gsize} --keep-dup=auto -n ${prefix}

    """
}
