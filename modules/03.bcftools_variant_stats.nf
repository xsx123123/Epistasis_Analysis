// 03.bcftools_variant_stats.nf 
process BCFTOOLS_VARIANT_STATS {
    tag "${sample_id}_BCFTOOLS_VARIANT_STATS"
    
    publishDir "${params.outdir}/03.variant_calling/variant_stats/${sample_id}", mode: 'link'
    
    input:
    tuple val(sample_id), path(vcf), path(csi), path(tbi)

    output:
    path "*.vcf.stats.txt", emit: variant_stats_report

    script:
    """
    bcftools stats \
        --threads ${task.cpus} \
        -f ${params.bcftools_reference} \
        ${vcf} > ${sample_id}.vcf.stats.txt
    """
}