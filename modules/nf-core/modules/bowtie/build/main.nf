process BOWTIE_BUILD {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f5ca09fd5aab931d9b87c532c69e0122ce5ff8ec88732f906e12108d48425e9/data' :
        'community.wave.seqera.io/library/bowtie_htslib_samtools:e1e242368ffcb5d3' }"

    input:
    path fasta

    output:
    path 'bowtie'       , emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p bowtie
    bowtie-build --threads ${task.cpus} ${fasta} bowtie/${fasta.baseName}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p bowtie
    touch bowtie/${fasta.baseName}.1.ebwt
    touch bowtie/${fasta.baseName}.2.ebwt
    touch bowtie/${fasta.baseName}.3.ebwt
    touch bowtie/${fasta.baseName}.4.ebwt
    touch bowtie/${fasta.baseName}.rev.1.ebwt
    touch bowtie/${fasta.baseName}.rev.2.ebwt
    """
}
