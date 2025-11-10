#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// param defaults can be found in nextflow.config
// module options can be found in conf/modules.config

// import modules
include { INPUT_CHECK } from './modules/local/input_check.nf'
include { generate_whitelist as whitelist_experimental ; generate_whitelist as whitelist_spikein } from './modules/local/generate_whitelist.nf'
include { FASTQC_TRIMGALORE } from './modules/local/fastqc_trimgalore.nf'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_EXPERIMENTAL ; BOWTIE2_BUILD as BOWTIE2_BUILD_SPIKEIN } from './modules/nf-core/modules/bowtie2/build/main.nf'
include { CHROMAP_INDEX as CHROMAP_INDEX_EXPERIMENTAL ; CHROMAP_INDEX as CHROMAP_INDEX_SPIKEIN } from './modules/nf-core/modules/chromap/index/main.nf'
include { CHROMAP_CHROMAP as CHROMAP_EXPERIMENTAL ; CHROMAP_CHROMAP as CHROMAP_SPIKEIN } from './modules/nf-core/modules/chromap/chromap/main.nf'
include { BOWTIE2_ALIGN as BOWTIE2_EXPERIMENTAL ; BOWTIE2_ALIGN as BOWTIE2_SPIKEIN } from './modules/nf-core/modules/bowtie2/align/main.nf'
include { samtools_whitelist_sort as samtools_experimental ; samtools_whitelist_sort as samtools_spikein } from './modules/local/samtools_whitelist_sort.nf'
include { PICARD_MERGESAMFILES as PICARD_MERGE_EXPERIMENTAL ; PICARD_MERGESAMFILES as PICARD_MERGE_SPIKEIN } from './modules/nf-core/modules/picard/mergesamfiles/main.nf'
include { PICARD_MARKDUPLICATES as PICARD_DEDUP_EXPERIMENTAL ; PICARD_MARKDUPLICATES as PICARD_DEDUP_SPIKEIN } from './modules/nf-core/modules/picard/markduplicates/main.nf'
include { overlap_check } from './modules/local/overlap_check.nf'
include { flagstat } from './modules/local/flagstat.nf'
include { calculate } from './modules/local/calculate.nf'
include { downsample } from './modules/local/downsample.nf'
include { macs2_peakcalling } from './modules/local/macs2_peakcalling.nf'
include { HOMER_ANNOTATEPEAKS } from './modules/nf-core/modules/homer/annotatepeaks/main.nf'
include { macs2_bdgcmp } from './modules/local/macs2_bdgcmp.nf'
include { MULTIQC } from './modules/local/multiqc.nf'
include { homer_findMotifsGenome } from './modules/local/homer_findMotifsGenome.nf'
include { idr } from './modules/local/idr.nf'
include { bedtools_consensus } from './modules/local/bedtools_consensus.nf'

workflow {

    // Check inputs and prepare FASTQ files for alignment:
    INPUT_CHECK (
        file(params.input),
        params.seq_center
    )

    FASTQC_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc,
        params.skip_trimming
    )

    // Assign experimental input files based on user parameters
    if (params.experimental == 'human') {
        ch_experimental_fa = params.human_fa
        ch_experimental_chro_index = params.human_chro_index ?: null
        ch_experimental_bowtie2_index = params.human_bowtie2_index ?: null
        ch_experimental_gtf = params.human_gtf ?: Channel.empty()
        ch_experimental_blacklist = params.human_blacklist ?: Channel.empty()
        ch_macs_gsize = params.macs_gsize_human
    }
    if (params.experimental == 'mouse') {
        ch_experimental_fa = params.mouse_fa
        ch_experimental_chro_index = params.mouse_chro_index ?: null
        ch_experimental_bowtie2_index = params.mouse_bowtie2_index ?: null
        ch_experimental_gtf = params.mouse_gtf ?: Channel.empty()
        ch_experimental_blacklist = params.mouse_blacklist ?: Channel.empty()
        ch_macs_gsize = params.macs_gsize_mouse
    }
    if (params.experimental == 'zebrafish') {
        ch_experimental_fa = params.zebrafish_fa
        ch_experimental_chro_index = params.zebrafish_chro_index ?: null
        ch_experimental_bowtie2_index = params.zebrafish_bowtie2_index ?: null
        ch_experimental_gtf = params.zebrafish_gtf ?: Channel.empty()
        ch_experimental_blacklist = params.zebrafish_blacklist ?: Channel.empty()
        ch_macs_gsize = params.macs_gsize_zebrafish
    }
    if (params.experimental == 'fly') {
        ch_experimental_fa = params.fly_fa
        ch_experimental_chro_index = params.fly_chro_index ?: null
        ch_experimental_bowtie2_index = params.fly_bowtie2_index ?: null
        ch_experimental_gtf = params.fly_gtf ?: Channel.empty()
        ch_experimental_blacklist = params.fly_blacklist ?: Channel.empty()
        ch_macs_gsize = params.macs_gsize_fly
    }

    // Assign spike-in input files based on user parameters
    if (params.spikein == 'human') {
        ch_spikein_fa = params.human_fa
        ch_spikein_chro_index = params.human_chro_index ?: null
        ch_spikein_bowtie2_index = params.human_bowtie2_index ?: null
        ch_spikein_blacklist = params.human_blacklist ?: Channel.empty()
    }
    if (params.spikein == 'mouse') {
        ch_spikein_fa = params.mouse_fa
        ch_spikein_chro_index = params.mouse_chro_index ?: null
        ch_spikein_bowtie2_index = params.mouse_bowtie2_index ?: null
        ch_spikein_blacklist = params.mouse_blacklist ?: Channel.empty()
    }
    if (params.spikein == 'zebrafish') {
        ch_spikein_fa = params.zebrafish_fa
        ch_spikein_chro_index = params.zebrafish_chro_index ?: null
        ch_spikein_bowtie2_index = params.zebrafish_bowtie2_index ?: null
        ch_spikein_blacklist = params.zebrafish_blacklist ?: Channel.empty()
    }
    if (params.spikein == 'fly') {
        ch_spikein_fa = params.fly_fa
        ch_spikein_chro_index = params.fly_chro_index ?: null
        ch_spikein_bowtie2_index = params.fly_bowtie2_index ?: null
        ch_spikein_blacklist = params.fly_blacklist ?: Channel.empty()
    }

    // Create indices for chosen aligner if not supplied
    if (params.aligner == 'chromap') {
        if (ch_experimental_chro_index == null) {
            CHROMAP_INDEX_EXPERIMENTAL (
                ch_experimental_fa
            )
            ch_experimental_chro_index = CHROMAP_INDEX_EXPERIMENTAL.out.index
        if (params.skip_downsample == false && ch_spikein_chro_index == null) {
            CHROMAP_INDEX_SPIKEIN (
                ch_spikein_fa
            )
            ch_spikein_chro_index = CHROMAP_INDEX_SPIKEIN.out.index
            }
        }
    }
    if (params.aligner == 'bowtie2') {
        if (ch_experimental_bowtie2_index == null) {
            BOWTIE2_BUILD_EXPERIMENTAL (
                ch_experimental_fa
            )
            ch_experimental_bowtie2_index = BOWTIE2_BUILD_EXPERIMENTAL.out.index
        }
        if (params.skip_downsample == false && ch_spikein_bowtie2_index == null) {
            BOWTIE2_BUILD_SPIKEIN (
                ch_spikein_fa
            )
            ch_spikein_bowtie2_index = BOWTIE2_BUILD_SPIKEIN.out.index
        }
    }
    
    // Create whitelist for samtools
    whitelist_experimental (
        ch_experimental_fa,
        ch_experimental_blacklist
    )

    whitelist_experimental
        .out
        .whitelist
        .set { ch_whitelist_experimental }

    // Align to experimental genome with either Chromap or Bowtie2
    if (params.aligner == 'chromap') {
            CHROMAP_EXPERIMENTAL (
                FASTQC_TRIMGALORE.out.reads,
                ch_experimental_fa,
                ch_experimental_chro_index,
                [], [], [], []
            )
    ch_aligned_experimental = CHROMAP_EXPERIMENTAL.out.bam
    }

    if (params.aligner == 'bowtie2') {
        BOWTIE2_EXPERIMENTAL (
            FASTQC_TRIMGALORE.out.reads,
            ch_experimental_bowtie2_index,
            [], []
        )
    ch_aligned_experimental = BOWTIE2_EXPERIMENTAL.out.bam
    }

    samtools_experimental (
        ch_aligned_experimental,
        ch_whitelist_experimental
    )

    // Group resequenced BAMs to prep for merge, adjust meta.id accordingly
    samtools_experimental
        .out
        .bam
        .map {
            meta, bam ->
                def meta_clone = meta.clone()
                meta_clone.remove('read_group')
                meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
                [ meta_clone, bam ] 
        }
        .groupTuple(by: [0])
        .map { 
            it ->
                [ it[0], it[1].flatten() ] 
        }
        .set { ch_sorted_experimental }

    PICARD_MERGE_EXPERIMENTAL (
        ch_sorted_experimental
    )

    PICARD_DEDUP_EXPERIMENTAL (
        PICARD_MERGE_EXPERIMENTAL.out.bam
    )

    // Spike-in specific steps (skipped if params.skip_downsample is true):
    if (params.skip_downsample != true) {

        whitelist_spikein (
            ch_spikein_fa,
            ch_spikein_blacklist
        )

        whitelist_spikein
            .out
            .whitelist
            .set { ch_whitelist_spikein }

        if (params.aligner == 'chromap') {
                CHROMAP_SPIKEIN (
                    FASTQC_TRIMGALORE.out.reads,
                    ch_spikein_fa,
                    ch_spikein_chro_index,
                    [], [], [], []
                )
        ch_aligned_spikein = CHROMAP_SPIKEIN.out.bam
        }

        if (params.aligner == 'bowtie2') {
            BOWTIE2_SPIKEIN (
                FASTQC_TRIMGALORE.out.reads,
                ch_spikein_bowtie2_index,
                [], []
            )
        ch_aligned_spikein = BOWTIE2_SPIKEIN.out.bam
        }

        samtools_spikein (
            ch_aligned_spikein,
            ch_whitelist_spikein
        )

        // Group resequenced spike-in BAMs to prep for merge, adjust meta.id accordingly
        samtools_spikein
            .out
            .bam
            .map {
                meta, bam ->
                    def meta_clone = meta.clone()
                    meta_clone.remove('read_group')
                    meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
                    [ meta_clone, bam ] 
            }
            .groupTuple(by: [0])
            .map { 
                it ->
                    [ it[0], it[1].flatten() ] 
            }
            .set { ch_sorted_spikein }

        PICARD_MERGE_SPIKEIN (
            ch_sorted_spikein
        )

        PICARD_DEDUP_SPIKEIN (
            PICARD_MERGE_SPIKEIN.out.bam
        )
    
        // Join experimental and spike-in BAMs in one channel and check for overlap
        PICARD_DEDUP_EXPERIMENTAL
            .out
            .bam
            .join ( PICARD_DEDUP_SPIKEIN.out.bam )
            .set { ch_experimental_spikein_bams }

        overlap_check (
            ch_experimental_spikein_bams,
            params.override_spikeinfail
        )
        
        // Add in spike-in status from overlap_check process to meta map
        overlap_check
            .out
            .experimental_bam
            .join(overlap_check.out.spikein_status)
            .map {
                meta, bam, spikein_status ->
                    [meta + [ spikein: spikein_status ], bam]
            }
            .set { ch_experimental_bam }

        overlap_check
            .out
            .spikein_bam
            .join(overlap_check.out.spikein_status)
            .map {
                meta, bam, spikein_status ->
                    [meta + [ spikein: spikein_status ], bam]
            }
            .set { ch_spikein_bam }

        // Collect overlap_check stats for output to user
        overlap_check
            .out
            .overlap_report
            .collectFile(keepHeader: true, skip: 1, storeDir: "${params.outdir}/${params.aligner}/overlap_check/")

        // Check if IP and control are both ok to downsample
        ch_experimental_bam
            .combine(ch_experimental_bam)
            .map {
                meta1, bam1, meta2, bam2 ->
                    meta1.control == meta2.id ? [ meta1, bam1, meta2, bam2 ] : null
            }
            .map {
                meta1, bam1, meta2, bam2 ->
                    def meta1_clone = meta1.clone()
                    def meta2_clone = meta2.clone()
                    meta1_clone.spikein = "no"
                    meta2_clone.spikein = "no"
                    meta1.spikein != meta2.spikein ? [ meta1_clone, bam1, meta2_clone, bam2 ] : [ meta1, bam1, meta2, bam2 ]
            }
            .multiMap {
                meta1, bam1, meta2, bam2 ->
                    ip: [ meta1, bam1 ]
                    control: [ meta2, bam2 ]
            }
            .set { ch_experimental_bam_split }

        ch_spikein_bam
            .combine(ch_spikein_bam)
            .map {
                meta1, bam1, meta2, bam2 ->
                    meta1.control == meta2.id ? [ meta1, bam1, meta2, bam2 ] : null
            }
            .map {
                meta1, bam1, meta2, bam2 ->
                    def meta1_clone = meta1.clone()
                    def meta2_clone = meta2.clone()
                    meta1_clone.spikein = "no"
                    meta2_clone.spikein = "no"
                    meta1.spikein != meta2.spikein ? [ meta1_clone, bam1, meta2_clone, bam2 ] : [ meta1, bam1, meta2, bam2 ]
            }
            .multiMap {
                meta1, bam1, meta2, bam2 ->
                    ip: [ meta1, bam1 ]
                    control: [ meta2, bam2 ]
            }
            .set { ch_spikein_bam_split }

        // Mix IP and control into one channel again, then split based on spikein status
        ch_experimental_bam_split
            .ip
            .mix(ch_experimental_bam_split.control)
            .branch { meta, bam ->
                yes: meta.spikein == "yes"
                no: meta.spikein == "no"
            }
            .set { ch_experimental_bam_spikestatus }

        ch_spikein_bam_split
            .ip
            .mix(ch_spikein_bam_split.control)
            .branch { meta, bam ->
                yes: meta.spikein == "yes"
                no: meta.spikein == "no"
            }
            .set { ch_spikein_bam_spikestatus }

        flagstat (
            ch_spikein_bam_spikestatus.yes
        )

        calculate (
            flagstat.out.spikein_count.collectFile()
        )

        // Prep and join bam files with associated scaling factors and info for downsampling
        calculate
            .out
            .scaled
            .splitCsv( header:true )
            .map { row -> 
                tuple( row.ID, row.spikein_reads_mapped, row.scaling_factor ) 
            }
            .set { scaling_ch }

        ch_experimental_bam_spikestatus
            .yes
            .map { meta, bam -> 
                tuple( bam.name, meta, bam ) 
            }
            .join ( scaling_ch )
            .map { ID, meta, bam, spikein_reads_mapped, scaling_factor -> 
                tuple( meta, scaling_factor, bam ) 
            }
            .set { joined }

        downsample (
            joined,
            params.seed
        )

        downsample
            .out
            .reads
            .set { ch_downsampled }

        // Create channel combining bam files of an IP with its control
        ch_downsampled
            .combine(ch_downsampled)
            .map {
                meta1, bam1, meta2, bam2 ->
                    meta1.control == meta2.id ? [ meta1, bam1, bam2 ] : null
            }
            .set { ch_downsampled_ip_control_bam }

        ch_experimental_bam_spikestatus
            .no
            .set { ch_deduped_experimental }
    }
    
    // If skipping normalization, get deduplicated experimental alignments ready and empty/placeholder channels created for downstream analysis
    // Otherwise, collect deduplication metrics for spike-in
    if (params.skip_downsample == true) {
        PICARD_DEDUP_EXPERIMENTAL
            .out
            .bam
            .set { ch_deduped_experimental }

        ch_downsampled_ip_control_bam = Channel.empty()
        ch_spikein_picard_metrics = Channel.empty()
    } else {
        PICARD_DEDUP_SPIKEIN
            .out
            .metrics
            .map { meta, metrics -> metrics }
            .collect()
            .set { ch_spikein_picard_metrics }
    }
    
    // Analysis steps downstream of normalization:
    // Collect outputs for MultiQC
    // Set up empty optional channels, in case they aren't created
    ch_FASTQC_metrics = Channel.empty()
    ch_Bowtie2_Experimental_metrics = Channel.empty()
    ch_Bowtie2_Spikein_metrics = Channel.empty()

    FASTQC_TRIMGALORE
        .out
        .fastqc_zip
        .map { meta, metrics -> metrics }
        .collect()
        .set { ch_FASTQC_metrics }

    PICARD_DEDUP_EXPERIMENTAL
        .out
        .metrics
        .map { meta, metrics -> metrics }
        .collect()
        .set { ch_experimental_picard_metrics }

    if (params.aligner == 'bowtie2') {
        BOWTIE2_EXPERIMENTAL
            .out
            .log
            .map { meta, metrics -> metrics }
            .collect()
            .set { ch_Bowtie2_Experimental_metrics }

        if (params.skip_downsample == false) {
            BOWTIE2_SPIKEIN
                .out
                .log
                .map { meta, metrics -> metrics }
                .collect()
                .set { ch_Bowtie2_Spikein_metrics }
        }
    }

    MULTIQC (
        ch_FASTQC_metrics.ifEmpty([]),
        ch_Bowtie2_Experimental_metrics.ifEmpty([]),
        ch_Bowtie2_Spikein_metrics.ifEmpty([]),
        ch_experimental_picard_metrics,
        ch_spikein_picard_metrics.ifEmpty([])
    )

    // Create channel combining bam files of an IP with its control for non-downsampled samples
    ch_deduped_experimental
        .combine(ch_deduped_experimental)
        .map {
            meta1, bam1, meta2, bam2 ->
                meta1.control == meta2.id ? [ meta1, bam1, bam2 ] : null
        }
        .set { ch_fullsize_ip_control_bam }

    // Mix downsampled and fullsize channels and send to downstream analysis modules
    ch_fullsize_ip_control_bam
        .mix(ch_downsampled_ip_control_bam)
        .set{ ch_ip_control_bam }

    // Call peaks using MACS2 subcommands
    macs2_peakcalling (
        ch_ip_control_bam,
        ch_macs_gsize,
        params.skip_downsample,
        params.macs2_cutoff,
        params.macs2_peak_method
    )

    if (params.skip_annotation == false) {
        HOMER_ANNOTATEPEAKS (
            macs2_peakcalling.out.peak,
            ch_experimental_fa,
            ch_experimental_gtf
        )
    }

    // Create bigWig files for visualization
    macs2_bdgcmp (
        macs2_peakcalling.out.chip_ctrl_bdg,
        whitelist_experimental.out.sizes,
        params.macs2_bigwig_method
    )

    if (params.skip_motif == false) {
        homer_findMotifsGenome (
            macs2_peakcalling.out.peak,
            ch_experimental_fa
        )
    }
    
    // Combine peaks based on input CSV's antibody column (in preparation for IDR)
    // IDR framework can only run on two peak files
    if (params.skip_idr == false && params.narrow_peak == true) {
        macs2_peakcalling
            .out
            .peak
            .filter { meta, peaks -> peaks.size() > 0 }
            .map { 
                    meta, peak -> 
                        [ meta.antibody, peak ] 
                }
            .groupTuple()
            .map {
                antibody, peaks ->
                    def meta_new = [:]
                    meta_new.id = antibody
                    [ meta_new, peaks ] 
            }
            .set { ch_combined_peaks }

        ch_combined_peaks
            .map { 
                meta, peaks ->
                    def count = peaks.size()
                    if (count != 2) {
                        log.warn "IDR can only be run on 2 peak files. Found ${count} for antibody/treatment \"${meta.id}\" - Skipping IDR for these samples."
                    }
                    count == 2 ? [ meta, peaks ] : null
            }
            .set { ch_combined_idr_peaks }
            
        idr (
            ch_combined_idr_peaks,
            params.idr_cutoff
        )
    }

    // Combine peaks based on input CSV's antibody column, then find simple consensus peak list
    // NOTE: Unlike IDR, can be run on more than two peak files
    if (params.skip_consensus == false) {
        macs2_peakcalling
            .out
            .peak
            .filter { meta, peaks -> peaks.size() > 0 }
            .map { 
                    meta, peak -> 
                        [ meta.antibody, peak ] 
                }
            .groupTuple()
            .map {
                antibody, peaks ->
                    def meta_new = [:]
                    meta_new.id = antibody
                    [ meta_new, peaks ] 
            }
            .set { ch_combined_peaks }

        bedtools_consensus (
            ch_combined_peaks
        )
    }
}

workflow.onComplete {
    println "\n"
    println "\033[1;34mPerCell pipeline completed at: $workflow.complete after $workflow.duration seconds.\033[0m"
    println "\033[1;34mExecution status: ${ workflow.success ? 'Completed successfully!' : 'Failed :/' }\033[0m"
}