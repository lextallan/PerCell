#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --job-name=RunPerCell

nextflow run PC.nf \
    -profile singularity \
    --human_fa '/path/to/downloaded/genome.fa' \
    --mouse_fa '/path/to/downloaded/genome.fa'