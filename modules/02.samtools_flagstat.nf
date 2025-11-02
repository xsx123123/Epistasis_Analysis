// 02.samtools_flagstat.nf
process SAMTOOLS_FLAGSTAT {
    tag "${sample_id}_FLAGSTAT"
    
    publishDir "${params.outdir}/02.mapping/${sample_id}", mode: 'link'
    
    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "*.flagstat.tsv", emit: flagstat_report

    script:
    """
    samtools flagstat \
             -@ ${task.cpus} \
             -O tsv \
             ${bam} > ${sample_id}.flagstat.tsv
    """
}
