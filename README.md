# PerCell ChIP-seq Pipeline
<br>
<img src="https://github.com/lextallan/PerCell/blob/main/subway.png" width="1250">
<br>

***PerCell ChIP-seq assays and their analysis using our pipeline are expected to generate similar qualitative readouts of protein:DNA interactions or histone modification localization to conventional ChIP-seq workflows, while additionally providing relative quantifications for these data across distinct replicates, lab groups, cellular treatment conditions, and genetic backgrounds.***


________________________________________

### Pipeline Setup and Execution

**<ins>Required Hardware:</ins>**

- Computer with Unix-based operating system

- A computer cluster or cloud-based computing platform for executing the pipeline (currently the pipeline is tested for compatibility with SLURM-managed platforms)

**<ins>Required Software:</ins>**

- Nextflow installation: [https://www.nextflow.io/docs/latest/install.html](https://www.nextflow.io/docs/latest/install.html)

- Singularity installation: [https://docs.sylabs.io/guides/latest/admin-guide/installation.html](https://docs.sylabs.io/guides/latest/admin-guide/installation.html)


**<ins>Data setup and automated pipeline execution</ins>**

*Timing: 30 min-1 h setup; 2-8 h execution*

1. Install both Nextflow and Singularity. Ensure both are properly running on the system.

2. If necessary, download files for chosen experimental and spike-in reference genomes in fasta format. Example commands for downloading and uncompressing human (hg38), mouse (mm10), zebrafish (danRer11), and fly (dm6) from UCSC can be found below:

    ```
    $ wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    $ gunzip hg38.fa.gz

    $ wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
    $ gunzip mm10.fa.gz

    $ wget http://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz
    $ gunzip danRer11.fa.gz`

    $ wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
    $ gunzip dm6.fa.gz
    ```

	
3. Create a csv file with names of samples, full locations of fastq files, antibodies/treatment used, and associated control sample names. Example input csv files are available in the examples directory. The file should be formatted as follows:

    - First line serves as column headings and must consist of the following: “sample,fastq_1,fastq_2,antibody,control”
    - Each additional line follows that format for each sequenced sample
    - Sample names should be unique, unless the same sample has been re-sequenced. Additional sequencing of the same sample can use identical sample names and the pipeline will merge the associated fastq files before alignment.
    - Samples with multiple replicates should have sample names ending in a similar pattern followed by the number of the replicate, e.g. “_1”/ “_2” or “_rep1”/ “_rep2”. These will be treated as replicates by the pipeline for downstream analyses.
    - fastq_1 and fastq_2 columns should contain the absolute path to the files, in an accessible directory. Files may be gzipped.
    - Antibody column should contain info about the antibody used as well as any treatments done on the sample. Samples with identical antibody columns will be combined for peak calling. For example, samples created using antibody A but distinct cell lines/ treatments B or C could be written as “antibodyA_B” and “antibodyA_C”, respectively. For input samples, this column should be left blank.
    - Control column should contain the exact name of the sample for the associated input, including replicate information. For example, “sampleA-Input_rep1” or “sampleA-Input_rep2” and NOT “sampleA-Input”.

4. OPTIONAL: Edit the file named `nextflow.config` with parameter choices (see execution step below for specific explanations) and options specific to the computing environment used (many organizations have suggested nextflow setup and config options, as collected by the nf-core community: https://nf-co.re/configs). The PerCell config file can be manually adjusted after directly downloading from this GitHub repository using the command

    `$ git clone https://raw.githubusercontent.com/lextallan/PerCell/master/nextflow.config`

    **Alternatively, parameter options can be specified at the command line at the point of pipeline execution.**

5. Execute the main pipeline script, PC.nf, with desired parameters (explained below). For ease of use, the entire pipeline and all potential parameters can be launched and customized from this single command. Make careful note of whether one or two dashes are used in the following parameters. Single dashes are for general nextflow options, double dashes are for pipeline-specific parameter choices.

    `$ nextflow run lextallan/PerCell/PC.nf`

    ```
    -profile singularity (Required. Tells nextflow to download/ use singularity images for each process)

    -qs: optional nextflow parameter to limit queue size, or number of processes that will be submitted at a time (e.g. to prevent overwhelming a shared cluster)

    --input: path to csv input file (see above for required format)

    --outdir: path to directory where the pipeline’s output will be found

    --aligner: choice of alignment tool with which samples will be aligned, <bowtie2,chromap> (default: bowtie2)

    --experimental: species identify from which the experimental cells originate, <human,mouse,zebrafish,fly> (default: human)

    --spikein: species identify from which the cells used for spike-in originate, <human,mouse,zebrafish,fly> (default: mouse)

    --human_fa: path to human reference fasta file, if used as either experimental or spike-in genome

    --mouse_fa: path to mouse reference fasta file, if used as either experimental or spike-in genome

    --zebrafish_fa: path to zebrafish reference fasta file, if used as if used as either experimental or spike-in genome

    --fly_fa: path to fly reference fasta file, if used as if used as either experimental or spike-in genome

    --skip_fastqc: choice of whether or not to skip FastQC quality control assessments, <true,false> (default: false)

    --skip_trimming: set to true if fastq files have already been trimmed of adapter sequences with an algorithm such as TrimGalore, <true,false> (default: false)

    --save_trimmed: if set to true, trimmed fastq files will be saved as output, <true,false> (default: false)

    --override_spikeinfail: the pipeline will not attempt to downsample experimental samples if spike-in content is detected to be below 0.5% of properly aligned total reads, setting this parameter to true will ignore this threshold and attempt normalization regardless of detected spike-in reads, <true,false> (default: false)

    --skip_downsample: set to true to disable PerCell normalization and carry out a standard ChIP-seq pipeline analysis, <true,false> (default: false)

    --seed: random seed used to downsample reads (default: 0)

    --macs2_peak_method: choice of method with which to score macs2 bedgraphs for peak calling, <ppois,qpois,subtract,logFE,FE,logLR,slogLR,max> (default: ppois)

    --macs2_cutoff: p/q-value cutoff used for calling peaks; unlike macs2’s default ‘callpeak’ command, cutoffs are taken in -log10(x) form. If using IDR, a particularly relaxed cutoff is highly recommended. For ChIPs targeting broader histone marks or when there are >2 replicates, IDR is not recommended and a more stringent cutoff of 1.301 (p/q-value of 0.05) or even 2 (p/q-value of 0.01) is suggested instead (default: 0.2218)

    --macs2_bigwig_method: choice of method with which to generate bigWig tracks for visualizing ChIP data within a genome browser, <ppois,qpois,subtract,logFE,FE,logLR,slogLR,max> (default: ppois)

    --skip_idr: if the IDR (Irreproducible Discovery Rate, https://github.com/nboley/idr) framework should be used to statistically identify significant peaks across two replicates. According to authors, the software is not designed for use with particularly broad peaks (e.g. those call in a ChIP for H3K9me3) and requires using a relaxed macs2 cutoff value in order to properly function, <true,false> (default: true)

    --idr_cutoff: cutoff value over which peaks will not be included within the output (default: 0.05)

    --skip_consensus: choice of whether consensus peaks should be identified across replicates, especially recommended along with a more stringent macs2 cutoff if there are >2 replicates, <true,false> (default: false)

    --skip_annotation: whether or not to use HOMER’s ‘annotatePeaks’ tool for annotation of peaks with nearest genes and genomic features <true,false> (default: false)

    --skip_motif: whether or not to use HOMER’s ‘findMotifsGenome’ tool to find enriched DNA sequence motifs in called peaks <true,false> (default: false)
    ```

________________________________________

### Summary of Pipeline Workflow

All pipeline outputs are stored in a directory named after the choice of alignment tool (`--aligner`) within the user-specified output directory (`--outidir`). When the applicable options are enabled (i.e. when `--skip_trimming` and `--skip_fastqc` set to false), adapter sequences are trimmed from reads and quality control data is collected and made available to the user via the MultiQC software in readable html format (MultiQC-Report-for-PerCell-pipeline.html). For user convenience, trimmed fastq files can be retained and organized for future use via the `--save_trimmed` parameter and stored in the trimgalore directory. 
Each fastq sequencing file is separately aligned to both the chosen experimental (`--experimental`) and spike-in (`--spikein`) genomes in a parallel fashion, with these intermediate alignments available in the aligned directory. Misaligned and duplicated reads are removed using the Samtools and Picard software suites, with post-filtering alignments and metrics data stored in the deduplicated directory. Summaries of these filtering metrics can be found in the picard_summaries directory. 

Next, the percentages of overlapping reads (independently aligning to both genomes) and spike-in reads (as a fraction of all aligned reads) are calculated for each sample, with the results collected in a single csv file (overlap_report.csv). Since the accuracy of normalization can be difficult when the percentage of spike-in reads is too low, the pipeline is designed to automatically exclude any samples with a calculated percentage of spike-in reads below 0.5%. This functionality can be overridden, forcing the pipeline to attempt normalization for all samples, by setting the `--override_spikeinfail` parameter to true. Conversely, PerCell normalization can be disabled (e.g. for analysis of ChIP-seq data lacking spike-in) by setting the parameter `--skip_downsample` to true. 

If PerCell normalization is enabled, scaling factors are calculated for each sample (scaling_factors.csv) based on the relative number of properly aligned spike-in reads. These factors are used to randomly (via the `--seed` parameter) downsample reads, with the now normalized alignment files stored in the downsampled directory. Each immunoprecipitated sample is matched with its corresponding control, and peaks are called using the MACS2 software. This is performed using MACS2’s subcommands to prevent the default MACS2 peak calling command’s normalization from obscuring our own. Each sample’s immunoprecipitated/control pairing is used to score potential peaks using the selected method (`macs2_peak_method`) and stored in the macs2_subcommand directory along with narrowPeak files containing peaks passing the cutoff given in the `--macs2_cutoff` parameter. Bedgraphs for each pair of immunoprecipitated and control samples are next compared using MACS2’s bdgcmp command based on the `--macs2_bigwig_method` parameter, before conversion into the visualizable bigWig format using tools from the UCSC software suite. 

Optionally, setting `--skip_idr` to false will use the IDR framework to identify statistically consistent peaks across two replicates (output in the idr directory) based on the chosen  `--idr_cutoff` value. As an alternative, e.g. when a more stringent MACS2 peak calling cutoff is used, consensus peaks can instead be identified using BEDTools by setting `--skip_consensus` to false. Finally, the HOMER software suite can be used to annotate peaksets via the `--skip_annotation` and identify enriched motifs via `--skip_motif` parameters, with the output found within their respective subdirectories inside the homer directory. 

________________________________________

### Output Directory Structure
```
<aligner>
    |-- aligned
    |   |-- <experimental_species>
    |   |   └── *.bam
    |   └── <spike-in_species>
    |       └── *.bam
    |-- deduplicated
    |   |-- <experimental_species>
    |   |   |-- *.dd.bam
    |   |   └── metrics
    |   |       └── *.metrics.txt
    |   └── <spike-in_species>
    |       |-- *.bam
    |       └── metrics
    |           └── *.metrics.txt
    |-- downsampled
    |   └── *_ds.bam
    |-- homer
    |   |-- annotations
    |   |   └── *.txt
    |   └── motifs
    |       └── *.html
    |-- macs2_bigwigs
    |   └── <macs2_bigwig_method>
    |       └── *.bigWig
    |-- macs2_subcommands
    |   |-- <macs2_cutoff>_cutoff
    |   |   └── *_subcommands.narrowPeak
    |   |-- <macs2_peak_method>_scored-bdg
    |   |   └── *.ppois-scored.bdg
    |   |-- <idr_cutoff>_idr
    |   |   └── *.narrowPeak
    |   └── consensus_peaks
    |       └── *.narrowPeak
    |-- multiqc
    |   └── MultiQC-Report-for-PerCell-pipeline.html
    |-- overlap_check
    |   └── overlap_report.csv
    |-- picard_summaries
    |   |-- <experimental_species>_summary.txt
    |   └── <spike-in_species>_summary.txt
    |-- scaling
        └── scaling_factors.csv
    └── trimgalore
        └── *.fq.gz
```

________________________________________

### Citations and Acknowledgements

A full list of references for software tools used within this pipeline can be found in the accompanying manuscript.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
 
> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.