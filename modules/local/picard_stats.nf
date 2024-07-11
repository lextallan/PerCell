process picard_stats {
    label 'process_medium'

    // custom Picard stat collection

    input:
    path experimental_picard_metrics, stageAs: 'experimental/*'
    path spikein_picard_metrics, stageAs: 'spikein/*'
    
    output:
    path "*_summary.txt", emit: stats

    script:
    """
    echo "sample,percent_duplication" > ${params.experimental}_summary.txt
    for file in ${experimental_picard_metrics}; do
        filename=\$(basename \$file)
        printf \$filename | awk -F. '{printf(\$1",")}' >> ${params.experimental}_summary.txt
        grep -A1 "^LIBRARY" \$file | tail -1 | cut -f 9 >> ${params.experimental}_summary.txt
    done

    echo "sample,percent_duplication" > ${params.spikein}_summary.txt
    for file in ${spikein_picard_metrics}; do
        filename=\$(basename \$file)
        printf \$filename | awk -F. '{printf(\$1",")}' >> ${params.spikein}_summary.txt
        grep -A1 "^LIBRARY" \$file | tail -1 | cut -f 9 >> ${params.spikein}_summary.txt
    done
    """
}
