process homer_findMotifsGenome {
    tag "$meta.id"
    label 'process_high'

    // Currently arguments in modules.config have this process set to only do known motifs.

    conda (params.enable_conda ? "bioconda::homer=4.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'quay.io/biocontainers/homer:4.11--pl526hc9558a2_3' }"

    input:
    tuple val(meta), path(peak)
    path  fasta

    output:
    path "knownResults.html", emit: known
    //path "output/homerResults.html", emit: denovo

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    findMotifsGenome.pl \\
        $peak \\
        $fasta \\
        . \\
        $args \\
        -p ${task.cpus}
    """
}
