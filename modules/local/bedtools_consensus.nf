process bedtools_consensus {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(combined_peaks)

    output:
    path '*.consensus.bed', emit: consensus_bed


    script:
    length = combined_peaks.size()
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${combined_peaks} | sort -k1,1 -k2,2n | bedtools merge -i stdin > merged_peaks.bed

    files=(*.narrowPeak)
    first=\${files[0]}
    bedtools intersect -wa -a merged_peaks.bed -b \${first} > running_consensus_\${first}.bed
    last=\${first}

    for i in "\${files[@]:1}"
    do
        bedtools intersect -wa -a running_consensus_\${last}.bed -b \$i > running_consensus_\${i}.bed
        last=\${i}
    done

    cp running_consensus_\${i}.bed ${prefix}.consensus.bed
    """
}
