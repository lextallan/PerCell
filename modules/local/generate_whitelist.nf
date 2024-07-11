process generate_whitelist {
    tag "$genome"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0':
        'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    path genome
    path blacklist

    output:
    path '*.bed', emit: whitelist

    script: // this custom module prepares a whitelist for use with filtering bam files
    """
    samtools faidx ${genome}
    cut -f1,2 ${genome}.fai > ${genome}.sizes
    sortBed -i ${blacklist} -g ${genome}.sizes | complementBed -i stdin -g ${genome}.sizes > ${genome.simpleName}_whitelist.bed
    """
}
