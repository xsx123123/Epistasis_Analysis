// 02.samtools_stats.nf
process SAMTOOLS_STATS {
    tag "${sample_id}_samtools_stats"
    
    publishDir "${params.outdir}/02.mapping/samtools_stats/${sample_id}", mode: 'link'
    
    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "*.stats.tsv", emit: stats_report

    script:
    """
    samtools stats \
             -@ ${task.cpus} \
             --reference ${params.bcftools_reference} \
             ${bam} > ${sample_id}.stats.tsv
    """
}