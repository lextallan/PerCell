process overlap_check {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(experimental_bam, stageAs: 'experimental/*'), path(spikein_bam, stageAs: 'spikein/*')
    val override

    output:
    tuple val(meta), path(experimental_bam), emit: experimental_bam
    tuple val(meta), path(spikein_bam), emit: spikein_bam
    tuple val (meta), env(spikein_used), emit: spikein_status
    path "overlap_report.csv", emit: overlap_report

    script:
    """
    samtools view $experimental_bam | cut -f 1 | LC_ALL=C sort | uniq > ${experimental_bam.baseName}.experimental.txt
    samtools view $spikein_bam | cut -f 1 | LC_ALL=C sort | uniq > ${spikein_bam.baseName}.spikein.txt
    LC_ALL=C comm -12 ${experimental_bam.baseName}.experimental.txt ${spikein_bam.baseName}.spikein.txt > ${meta.id}_common.txt

    common=\$(cat ${meta.id}_common.txt | wc -l)
    experimental=\$(cat ${experimental_bam.baseName}.experimental.txt | wc -l)
    spikein=\$(cat ${spikein_bam.baseName}.spikein.txt | wc -l)

    echo "Common reads: \$common"
    echo "Experimental reads: \$experimental"
    echo "Spikein reads: \$spikein"

    percent_spiked=\$(awk -v spikein="\$spikein" -v experimental="\$experimental" 'BEGIN { OFMT="%.4f"; x= spikein / experimental; print x}')
    echo "Percent spiked-in: \${percent_spiked}"

    percent_common=\$(awk -v common="\$common" -v experimental="\$experimental" 'BEGIN { OFMT="%.4f"; x= common / experimental; print x}')
    echo "Percent_common: \${percent_common}"
    
    # if the percent spiked in experimental and spike-in samples is greater than 0.5% then we consider a spike-in as used; 0 or 1 values passed through channel, then summed to check if every sample has spike-in content
    if [[ \${percent_spiked} > 0.02 ]]
    then
        spikein_used="yes"
    elif  [ $override == true ]
    then
        spikein_used="yes"
    else
        spikein_used="no"
    fi

    echo "Spikein_used: \${spikein_used}"

    if [ $override == true ]
    then
        echo "WARNING: Overriding spike-in status, bams will be downsampled even if very few spike-in reads detected!"
    fi

    # Collect data for each sample -- write to csv and collect to report to user
    echo "Sample,Experimental_reads,Spikein_reads,Percent_spiked_in,Spikein_used,Common_reads,Percent_common" > overlap_report.csv
    echo "${meta.id},\$experimental,\$spikein,\${percent_spiked},\${spikein_used},\$common,\${percent_common}" >> overlap_report.csv

    """
}
