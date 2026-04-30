process deeptools_bamCoverage {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::deeptools=3.5.1 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:2c687053c0252667cca265c9f4118f2c205a604c-0':
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:2c687053c0252667cca265c9f4118f2c205a604c-0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig  , optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"

    """
    samtools index $input
    
    bamCoverage \\
        --bam $input \\
        $args \\
        --numberOfProcessors ${task.cpus} \\
        --outFileName ${prefix}.bigWig
    """

}