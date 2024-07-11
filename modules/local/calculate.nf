process calculate {
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
         'https://depot.galaxyproject.org/singularity/mulled-v2-6f885f082106c9deee6bb625a3b516488915a563:f51a296ed289dee746efb7f2fdd520f645844f85-0' :
         'quay.io/biocontainers/mulled-v2-6f885f082106c9deee6bb625a3b516488915a563:f51a296ed289dee746efb7f2fdd520f645844f85-0' }"

    input:
    path spikein_count

    output:
    path 'scaling_factors.csv', emit: scaled

    """
    calculate_min_spikein.py ${spikein_count}
    """
}
