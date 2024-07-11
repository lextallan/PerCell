process downsample {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), val(scaling_factor), path(bam)
    val seed

    output:
    tuple val(meta), path("*_ds.bam"), emit: reads

    script:
    def seed_sample = seed ? "--subsample-seed $seed" : ''
    """
    samtools view -h -b -o ${meta.id}_ds.bam -f 1 --subsample ${scaling_factor} ${seed_sample} ${bam}
    """
}
