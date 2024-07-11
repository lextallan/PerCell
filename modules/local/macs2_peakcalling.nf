process macs2_peakcalling {
    tag "$meta.id"
    label 'process_high_memory'

    conda (params.enable_conda ? "bioconda::macs2=2.2.9.1, bioconda::mawk=1.3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-447fd94c2e07f56860fcc889687065479b8d5c34:c8f8a83fdae750b7e437ecfac55a5ef14b97fc2b-0' :
        'quay.io/biocontainers/mulled-v2-447fd94c2e07f56860fcc889687065479b8d5c34:c8f8a83fdae750b7e437ecfac55a5ef14b97fc2b-0' }"

    input:
    tuple val(meta), path(ipbam), path(controlbam)
    val macs2_gsize
    val skip_downsample
    val macs2_cutoff
    val macs2_method

    output:
    tuple val(meta), path("*_subcommands.narrowPeak"), emit: peak
    tuple val(meta), path("chip_bedgraph.bdg"), path("control_bedgraph.bdg"), emit: chip_ctrl_bdg
    tuple val(meta), path("*-scored.bdg"), emit: scored_bdg
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    macs2 filterdup -i $ipbam --keep-dup all -o ${ipbam.baseName}.bed
    macs2 filterdup -i $controlbam --keep-dup all -o ${controlbam.baseName}.bed

    # USE MACS2's callpeak to get fragment length d and half fragment length d2
    macs2 callpeak --treatment $ipbam --control $controlbam -f BAMPE -g $macs2_gsize -n ${prefix}_callpeak
    d=\$(cat ${prefix}_callpeak_peaks.xls | awk '/# d/ {print \$4}')
    d2=\$(awk -v d="\$d" 'BEGIN { x= d / 2; print int(x)}')
    echo "Fragment length: \$d"
    echo "Half fragment length: \$d2"

    #Convert control bed file into single ended bed (READ LENGTH set to 150)
    mawk -v OFS="\\t" -v RL=150 '{print \$1,\$2,\$2+RL,".",".","+\\n"\$1,\$3-RL,\$3,".",".","-"}' \\
        ${controlbam.baseName}.bed > ${controlbam.baseName}_as_SE.bed

    macs2 \\
        pileup \\
        -f BEDPE \\
        -i ${ipbam.baseName}.bed \\
        -o ${ipbam.baseName}.pileup.bdg

    # d bkgrd
    macs2 \\
        pileup \\
        -f BED \\
        -i ${controlbam.baseName}_as_SE.bed \\
        -B \\
        -o ${controlbam.baseName}.d_bkgrd.bdg \\
        --extsize \$d2

    # small local bkgrd (1k, extsize 500)
    macs2 \\
        pileup \\
        -f BED \\
        -i ${controlbam.baseName}_as_SE.bed \\
        -B \\
        -o ${controlbam.baseName}.s_bkgrd.bdg \\
        --extsize 500 

    # large local bkgrd (10k, extsize 5000)
    macs2 \\
        pileup \\
        -f BED \\
        -i ${controlbam.baseName}_as_SE.bed \\
        -B \\
        -o ${controlbam.baseName}.l_bkgrd.bdg \\
        --extsize 5000

    # Halve lambda values 
    mawk -v OFS="\\t" '{print \$1,\$2,\$3,\$4/2}' ${controlbam.baseName}.d_bkgrd.bdg > ${controlbam.baseName}.d_bkgrd.halved.bdg
    mawk -v OFS="\\t" '{print \$1,\$2,\$3,\$4/2}' ${controlbam.baseName}.s_bkgrd.bdg > ${controlbam.baseName}.s_bkgrd.halved.bdg
    mawk -v OFS="\\t" '{print \$1,\$2,\$3,\$4/2}' ${controlbam.baseName}.l_bkgrd.bdg > ${controlbam.baseName}.l_bkgrd.halved.bdg

    # Normalize slocal
    ps=\$(awk -v d="\$d" 'BEGIN { OFMT="%.4f"; x= d / 1000; print x}')
    macs2 bdgopt -i ${controlbam.baseName}.s_bkgrd.halved.bdg -m multiply -p \$ps -o ${controlbam.baseName}.s_bkgrd.halved.norm.bdg

    # Normalize llocal
    pl=\$(awk -v d="\$d" 'BEGIN { OFMT="%.4f"; x= d / 10000; print x}')
    macs2 bdgopt -i ${controlbam.baseName}.l_bkgrd.halved.bdg -m multiply -p \$pl -o ${controlbam.baseName}.l_bkgrd.halved.norm.bdg

    # Combine s and l local max
    macs2 bdgcmp -m max -t ${controlbam.baseName}.s_bkgrd.halved.norm.bdg -c ${controlbam.baseName}.l_bkgrd.halved.norm.bdg -o ${controlbam.baseName}.sl_bkgrd.halved.norm.bdg

    # Combine sl with d max
    macs2 bdgcmp -m max -t ${controlbam.baseName}.sl_bkgrd.halved.norm.bdg -c ${controlbam.baseName}.d_bkgrd.halved.bdg -o ${controlbam.baseName}.dsl_bkgrd.halved.norm.bdg

    # Genome wide bkgrd and combine
    ctrl_num_read=\$(wc -l ${controlbam.baseName}_as_SE.bed | awk '{print \$1}')
    p_genome=\$(awk -v d="\$d" -v g="$macs2_gsize" -v c="\$ctrl_num_read" 'BEGIN { OFMT="%.4f"; z= c * d; x= z / g; print x}')
    echo "p_genome: \$p_genome"
    macs2 bdgopt -i ${controlbam.baseName}.dsl_bkgrd.halved.norm.bdg -m max -p \$p_genome -o ${controlbam.baseName}.local_bias_raw.bdg

    # Scale treatment and control only if no downsampling was performed
    if [ "${meta.spikein}" == "no" ] || [ ${skip_downsample} == true ]
    then
        echo "No downsampling, so scaling treatment/control here instead"
        treat_num_read=\$(wc -l ${ipbam.baseName}.bed | awk '{print \$1}')
        ctrl_num_read=\$(wc -l ${controlbam.baseName}.bed | awk '{print \$1}')
        ratio=\$(awk -v t="\$treat_num_read" -v c="\$ctrl_num_read" 'BEGIN { OFMT="%.4f"; x= t / c; print x}')
        if [[ \${ratio} < 1 ]]
        then
            echo "Scaling control"
            macs2 bdgopt -i ${controlbam.baseName}.local_bias_raw.bdg -m multiply -p \$ratio -o ${controlbam.baseName}.local_bias.bdg
            macs2 bdgcmp -t ${ipbam.baseName}.pileup.bdg -c ${controlbam.baseName}.local_bias.bdg -m ${macs2_method} -o ${ipbam.baseName}.${macs2_method}-scored.bdg
            mv ${ipbam.baseName}.pileup.bdg chip_bedgraph.bdg
            mv ${controlbam.baseName}.local_bias_raw.bdg control_bedgraph.bdg
        else
            echo "Scaling treatment"
            inverse_ratio=\$(awk -v r="\$ratio" 'BEGIN { OFMT="%.4f"; x= 1 / r; print x}')
            macs2 bdgopt -i ${ipbam.baseName}.pileup.bdg -m multiply -p \${inverse_ratio} -o ${ipbam.baseName}.pileup_scaled.bdg
            macs2 bdgcmp -t ${ipbam.baseName}.pileup_scaled.bdg -c ${controlbam.baseName}.local_bias_raw.bdg -m ${macs2_method} -o ${ipbam.baseName}.${macs2_method}-scored.bdg
            mv ${ipbam.baseName}.pileup.bdg chip_bedgraph.bdg
            mv ${controlbam.baseName}.local_bias_raw.bdg control_bedgraph.bdg
        fi
    else
        echo "Samples were normalized via downsampling, so no scaling needed"
        macs2 bdgcmp -t ${ipbam.baseName}.pileup.bdg -c ${controlbam.baseName}.local_bias_raw.bdg -m ${macs2_method} -o ${ipbam.baseName}.${macs2_method}-scored.bdg
        mv ${ipbam.baseName}.pileup.bdg chip_bedgraph.bdg
        mv ${controlbam.baseName}.local_bias_raw.bdg control_bedgraph.bdg
    fi

    # NOTE: MACS2 callpeak's default cutoff value is 0.05, but bdgpeakcall takes it as -log10(cutoff), so -log10(0.05) = -c 1.301, -log10(0.01) = -c 2. Lenient cutoff for IDR: -log10(0.6) = -c 0.2218
    macs2 bdgpeakcall -i ${ipbam.baseName}.${macs2_method}-scored.bdg -c ${macs2_cutoff} -l \$d -g 150 -o ${ipbam.baseName}.peaks.bed
    tail +2 ${ipbam.baseName}.peaks.bed > ${ipbam.baseName}_subcommands.narrowPeak

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs2: \$(macs2 --version | sed -e "s/macs2 //g")
    END_VERSIONS
    """
}
