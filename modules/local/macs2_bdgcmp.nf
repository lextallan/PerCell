process macs2_bdgcmp {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::macs2=2.2.9.1 bioconda::bedtools=2.31.0 bioconda::ucsc-bedclip=377 bioconda::ucsc-bedgraphtobigwig=445" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-20f97261dc026feb7aca77ec7eca9ebfcb93f1ef:f4f30a4635214a5ded57a1b274e5f847eff9aa0b-0' :
        'quay.io/biocontainers/mulled-v2-20f97261dc026feb7aca77ec7eca9ebfcb93f1ef:f4f30a4635214a5ded57a1b274e5f847eff9aa0b-0' }"

    input:
    tuple val(meta), path(treatment_bdg), path(control_bdg)
    path chr_sizes
    val method

    output:
    tuple val(meta), path("*_signal.bigWig"), emit: bigwig
    path  "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def full_method = method == 'logFE' ? 'logFE -p 0.00001' : method == 'logLR' ? 'logLR -p 0.00001' : method
    def VERSION_BEDCLIP = '377' // WARN: Version information not available from ucsc CLI, update manually
    def VERSION_BEDGRAPHTOBIGWIG = '445' // WARN: Version information not available from ucsc CLI, update manually
    """
    macs2 \\
        bdgcmp \\
        -t $treatment_bdg \\
        -c $control_bdg \\
        -o ${prefix}_${method}.bdg \\
        -m ${full_method}

    LC_COLLATE=C sort -k1,1 -k2,2n ${prefix}_${method}.bdg > ${prefix}_${method}.sorted.bdg

    slopBed -i ${prefix}_${method}.sorted.bdg -g ${chr_sizes} -b 0 | \\
        awk '{if (\$3 != -1) print \$0}' | \\
        bedClip stdin ${chr_sizes} ${prefix}_signal.bedgraph

    # Merge bed to enable bedgraphtobigwig without error
    # !Determine if still necessary!
    mergeBed -i ${prefix}_signal.bedgraph -d -1 -c 4 -o max > ${prefix}_signal.merged.bedgraph

    bedGraphToBigWig ${prefix}_signal.merged.bedgraph ${chr_sizes} ${prefix}_signal.bigWig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs2: \$(macs2 --version | sed -e "s/macs2 //g")
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        bedClip: $VERSION_BEDCLIP
        bedGraphToBigWig: $VERSION_BEDGRAPHTOBIGWIG
    END_VERSIONS
    """
}
