process flagstat {
    tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    path 'spikein_mapped.csv', emit: spikein_count
    path '*.tsv', emit: stats

    """
    samtools flagstat ${bam} -O tsv > ${meta.id}.tsv
    grep 'primary mapped' ${meta.id}.tsv | grep -v '%'| cut -f1 | xargs echo "${bam}," > spikein_mapped.csv
    """
}
