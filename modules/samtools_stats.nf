process SAMTOOLS_STATS {
    tag "$sample_id"
    publishDir "${params.outdir}/03_mapping/samtools_stats/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "*.stats.tsv", emit: stats_report

    conda "${baseDir}/envs/bwa2.yaml" // samtools is in bwa2.yaml

    script:
    """
    samtools stats \
             -@ ${task.cpus} \
             --reference ${params.bcftools_reference} \
             ${bam} > ${sample_id}.stats.tsv
    """
}