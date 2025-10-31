process BCFTOOLS_CALL {
    tag "$sample_id"
    publishDir "${params.outdir}/04_variant_calling/${sample_id}", mode: 'copy', pattern: "*.raw.vcf.gz"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*.raw.vcf.gz"), emit: raw_vcf

    script:
    """
    ${params.bcftools_path} mpileup \
        --threads ${task.cpus} \
        --verbosity ${params.bcftools_verbosity} \
        -f ${params.bcftools_reference} \
        ${bam} | ${params.bcftools_path} call \
        --threads ${task.cpus} \
        --ploidy ${params.bcftools_ploidy} \
        --verbosity ${params.bcftools_verbosity} \
        -mv -Oz -o ${sample_id}.raw.vcf.gz
    """
}