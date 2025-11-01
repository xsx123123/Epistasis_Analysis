#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    Epistasis Analysis Pipeline
========================================================================================
*/

// Parameters are defined in nextflow.config

// ---------------- //
//    Help Message
// ---------------- //

if (params.help) {
    log.info '''
    Epistasis Analysis Pipeline
    =================================
    Usage:
    nextflow run main.nf -profile <docker/conda>

    Parameters:
      --sample_csv      Path to the sample sheet file (default: 'sample.csv')
      --raw_data_path   Path to the raw data directory
      --outdir          Path to the output directory (default: './results')
      --fasta           Path to the reference genome FASTA file (REQUIRED)
      --help            Display this help message
    '''.stripIndent()
    exit 0
}

// ---------------- //
//    Modules
// ---------------- //

// include nextflow modules
// raw data quality control & clean process
include { FASTQC_RAW } from './modules/01.fastqc_raw.nf'
include { FASTP_CLEAN } from './modules/01.fastp_clean.nf'
include { FASTQ_SCREEN } from './modules/01.fastq_screen.nf'
include { MULTIQC } from './modules/01.multiqc.nf'
// mapping process with bwa algorith
include { BWA_INDEX } from './modules/02.bwa_index.nf'
include { BWA_MEM } from './modules/02.bwa_mem.nf'
include { SAMTOOLS_SORT_INDEX } from './modules/02.samtools_sort_index.nf'
include { SAMBAMBA_MARKDUPLICATES } from './modules/02.sambamba_markduplicates.nf'
include { MULTIQC_MAPPING } from  './modules/02.multiqc_mapping.nf'
include { BAM_COVERAGE } from './modules/02.bam_coverage.nf'
include { QUALIMAP_QC } from './modules/02.qualimap_qc.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/02.samtools_flagstat.nf'
include { SAMTOOLS_STATS } from './modules/02.samtools_stats.nf'
// variant calling process
include { BCFTOOLS_CALL } from './modules/03.bcftools_call.nf'
include { BCFTOOLS_SORT_INDEX } from './modules/03.bcftools_sort_index.nf'
include { BCFTOOLS_MERGE } from './modules/03.bcftools_merge.nf'
include { BCFTOOLS_FILTER } from './modules/03.bcftools_filter.nf'
include { SNPEFF_ANNOTATION } from './modules/03.snpEff_annotation.nf'
include { BCFTOOLS_VARIANT_STATS } from './modules/03.bcftools_variant_stats.nf'
include { MULTIQC_VARIANT_STATS } from './modules/03.multiqc_variant_stats.nf'

// ---------------- //
//    Workflow
// ---------------- //

workflow {
    // Analysis init input channel
    ch_samples = Channel
    .fromPath(params.sample_csv) // loading nextflow channel from csv file
    .splitCsv(header:true)
    .map { row -> 
        tuple(
            row.sample_id, 
            file("${params.raw_data_path}/${row.fq1}"), 
            file("${params.raw_data_path}/${row.fq2}")
        )
    } // .view { "Sample tuple: $it" }
    
    // Quality control was performed on raw data using FastQC
    FASTQC_RAW(ch_samples)

    // FASTQ_SCREEN
    FASTQ_SCREEN(ch_samples)

    // Cleaning of raw data using FastP
    FASTP_CLEAN(ch_samples)

    // define raw & clean data channels
    def all_reports_for_multiqc = FASTQC_RAW.out.zip
                                .mix(FASTP_CLEAN.out.json)
                                .mix(FASTP_CLEAN.out.html)
                                .mix(FASTQ_SCREEN.out.screen_txt)
                                .mix(FASTQ_SCREEN.out.screen_html)

    def multiqc_input_list = all_reports_for_multiqc.collect()
    
    // MultiQC report for raw and clean data
    MULTIQC(multiqc_input_list)

    // define reference genome channel
    ch_fasta = Channel.fromPath(params.genome)
    
    // define bwa2-mem index & mapping channels
    ch_bwa_index = Channel.empty()
    if (params.bwa_index_dir) {
        ch_bwa_index = Channel.fromPath(params.bwa_index_dir)
        log.info "Using existing BWA index from: ${params.bwa_index_dir}"
    } else {
        BWA_INDEX(ch_fasta)
        ch_bwa_index = BWA_INDEX.out.index
    }
    
    // clean data mapping by bwa-mem2
    BWA_MEM(FASTP_CLEAN.out.trimmed_reads,
            ch_bwa_index)

    // sort and index the mapped bam file
    SAMTOOLS_SORT_INDEX(BWA_MEM.out.bam)

    // mark duplicates by samblba
    SAMBAMBA_MARKDUPLICATES(SAMTOOLS_SORT_INDEX.out.sorted_indexed_bam)

    // mapping result status channels
    BAM_COVERAGE(SAMBAMBA_MARKDUPLICATES.out.marked_duplicates_bam)

    // mapping quality control by QUALIMAP
    QUALIMAP_QC(SAMBAMBA_MARKDUPLICATES.out.marked_duplicates_bam)

    // mapping FLAGSTAT
    SAMTOOLS_FLAGSTAT(SAMBAMBA_MARKDUPLICATES.out.marked_duplicates_bam)

    // mapping STATS
    SAMTOOLS_STATS(SAMBAMBA_MARKDUPLICATES.out.marked_duplicates_bam)

    // define mapping stats channels
    def mapping_reports_for_multiqc = SAMTOOLS_STATS.out.stats_report
                                .mix(SAMTOOLS_FLAGSTAT.out.flagstat_report)
                                .mix(BAM_COVERAGE.out.global_dist)
                                .mix(BAM_COVERAGE.out.summary)

    def mapping_multiqc_input_list = mapping_reports_for_multiqc.collect()

    // MultiQC report for mapping result
    MULTIQC_MAPPING(mapping_multiqc_input_list)
    
    // Call variant
    BCFTOOLS_CALL(SAMBAMBA_MARKDUPLICATES.out.marked_duplicates_bam)
    
    // VCF Sort & Index
    BCFTOOLS_SORT_INDEX(BCFTOOLS_CALL.out.raw_vcf)

    // VCF Merge
    BCFTOOLS_MERGE(BCFTOOLS_SORT_INDEX.out.sorted_indexed_vcf.map{ it[1] }.collect(),
                   BCFTOOLS_SORT_INDEX.out.sorted_indexed_vcf.map{ it[2] }.collect(),
                   BCFTOOLS_SORT_INDEX.out.sorted_indexed_vcf.map{ it[3] }.collect())

    // VCF Filter
    BCFTOOLS_FILTER(
        BCFTOOLS_MERGE.out.merged_vcf,
        BCFTOOLS_MERGE.out.merged_vcf_csi,
        BCFTOOLS_MERGE.out.merged_vcf_tbi
    )

    // VCF Annotate
    SNPEFF_ANNOTATION(
        BCFTOOLS_FILTER.out.filtered_vcf,
        BCFTOOLS_FILTER.out.filtered_vcf_csi,
        BCFTOOLS_FILTER.out.filtered_vcf_tbi
    )

    // Variant stats by bcftools
    BCFTOOLS_VARIANT_STATS(BCFTOOLS_SORT_INDEX.out.sorted_indexed_vcf)

    // merge Variant stats report by multiqc
    MULTIQC_VARIANT_STATS(BCFTOOLS_VARIANT_STATS.out.variant_stats_report.collect())
}