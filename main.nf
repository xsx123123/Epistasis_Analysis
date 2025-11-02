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
      --help            Display this help message
    '''.stripIndent()
    exit 0
}

// ---------------- //
//    Modules
// ---------------- //

include { FASTQC_RAW } from './modules/fastqc_raw.nf'
include { FASTP_CLEAN } from './modules/fastp_clean.nf'
include { BWA_INDEX } from './modules/bwa_index.nf'
include { BWA_MEM } from './modules/bwa_mem.nf'
include { SAMTOOLS_SORT_INDEX } from './modules/samtools_sort_index.nf'
include { SAMBAMBA_MARKDUPLICATES } from './modules/sambamba_markduplicates.nf'
include { SAMTOOLS_INDEX_DUP } from './modules/samtools_index_dup.nf'
include { BCFTOOLS_CALL } from './modules/bcftools_call.nf'
include { BCFTOOLS_SORT_INDEX } from './modules/bcftools_sort_index.nf'
include { BCFTOOLS_MERGE } from './modules/bcftools_merge.nf'
include { BCFTOOLS_FILTER } from './modules/bcftools_filter.nf'
include { SNPEFF_ANNOTATION } from './modules/snpEff_annotation.nf'
include { BAM_COVERAGE } from './modules/bam_coverage.nf'
include { QUALIMAP_QC } from './modules/qualimap_qc.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/samtools_flagstat.nf'
include { SAMTOOLS_STATS } from './modules/samtools_stats.nf'
include { BCFTOOLS_VARIANT_STATS } from './modules/bcftools_variant_stats.nf'
include { MULTIQC_VARIANT_STATS } from './modules/multiqc_variant_stats.nf'

// ---------------- //
//    Workflow
// ---------------- //

workflow {
    ch_samples = Channel
        .fromPath(params.sample_csv)
        .splitCsv(header:true)
        .map { row -> tuple(row.sample_id, file(row.fq1), file(row.fq2)) }

    FASTQC_RAW(ch_samples)

    ch_trimmed_reads = FASTP_CLEAN(ch_samples).trimmed_reads

    ch_bwa_index = BWA_INDEX(file(params.genome))

    ch_mapped_reads = BWA_MEM(ch_trimmed_reads, ch_bwa_index).bam

    ch_sorted_indexed_bam = SAMTOOLS_SORT_INDEX(ch_mapped_reads).sorted_indexed_bam

    ch_marked_duplicates_bam = SAMBAMBA_MARKDUPLICATES(ch_sorted_indexed_bam).marked_duplicates_bam

    ch_indexed_dup_bam = SAMTOOLS_INDEX_DUP(ch_marked_duplicates_bam).indexed_dup_bam

    ch_raw_vcf = BCFTOOLS_CALL(ch_indexed_dup_bam).raw_vcf

    ch_sorted_indexed_vcf = BCFTOOLS_SORT_INDEX(ch_raw_vcf).sorted_indexed_vcf

    ch_merged_vcf_outputs = BCFTOOLS_MERGE(ch_sorted_indexed_vcf.collect())

    ch_merged_vcf = ch_merged_vcf_outputs.merged_vcf

    ch_merged_vcf_csi = ch_merged_vcf_outputs.merged_vcf_csi
    
    ch_merged_vcf_tbi = ch_merged_vcf_outputs.merged_vcf_tbi

    ch_filtered_vcf_outputs = BCFTOOLS_FILTER(ch_merged_vcf, ch_merged_vcf_csi, ch_merged_vcf_tbi)
    ch_filtered_vcf = ch_filtered_vcf_outputs.filtered_vcf
    ch_filtered_vcf_csi = ch_filtered_vcf_outputs.filtered_vcf_csi
    ch_filtered_vcf_tbi = ch_filtered_vcf_outputs.filtered_vcf_tbi

    ch_annotated_vcf_outputs = SNPEFF_ANNOTATION(ch_filtered_vcf, ch_filtered_vcf_csi, ch_filtered_vcf_tbi)
    ch_annotated_csv = ch_annotated_vcf_outputs.annotated_csv
    ch_annotated_html = ch_annotated_vcf_outputs.annotated_html
    ch_annotated_vcf = ch_annotated_vcf_outputs.annotated_vcf

    ch_coverage_outputs = BAM_COVERAGE(ch_indexed_dup_bam)
    ch_coverage_global_dist = ch_coverage_outputs.global_dist
    ch_coverage_summary = ch_coverage_outputs.summary

    ch_qualimap_outputs = QUALIMAP_QC(ch_indexed_dup_bam)
    ch_qualimap_html = ch_qualimap_outputs.html_report
    ch_qualimap_txt = ch_qualimap_outputs.txt_report

    ch_flagstat_report = SAMTOOLS_FLAGSTAT(ch_indexed_dup_bam).flagstat_report

    ch_samtools_stats_report = SAMTOOLS_STATS(ch_indexed_dup_bam).stats_report

    ch_variant_stats_report = BCFTOOLS_VARIANT_STATS(ch_sorted_indexed_vcf).variant_stats_report

    ch_multiqc_variant_stats_report = MULTIQC_VARIANT_STATS(ch_variant_stats_report.collect()).multiqc_report

    ch_multiqc_variant_stats_report.view()
}