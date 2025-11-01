// 03.bcftools_call.nf
process BCFTOOLS_CALL {
    
    tag "${sample_id}_Bcftools_call_variant"

    publishDir "${params.outdir}/03.variant_calling/${sample_id}", mode: 'link', pattern: "*.raw.vcf.gz"
    
    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*.raw.vcf.gz"), emit: raw_vcf

    script:
    """
    bcftools mpileup \
        --threads ${task.cpus} \
        --verbosity ${params.bcftools_verbosity} \
        -f ${params.bcftools_reference} \
        ${bam} | bcftools call \
        --threads ${task.cpus} \
        --ploidy ${params.bcftools_ploidy} \
        --verbosity ${params.bcftools_verbosity} \
        -mv -Oz -o ${sample_id}.raw.vcf.gz
    """
}