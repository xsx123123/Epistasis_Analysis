// 03.bcftools_call.nf
process BCFTOOLS_CALL_BY_CHR {
    tag "${sample_id}_${chr}_Bcftools_call"

    publishDir "${params.outdir}/03.variant_calling/${sample_id}_chr/", mode: 'link', pattern: "${sample_id}.${chr}.raw.vcf.gz"

    input:
    tuple val(sample_id), path(bam), path(bai), val(chr)  // 改为单个元组输入

    output:
    tuple val(sample_id), val(chr), path("*.raw.vcf.gz"), emit: vcf_by_chr

    script:
    def out_vcf = "${sample_id}.${chr}.raw.vcf.gz"
    def small_cpus = 2 
    
    """
    bcftools mpileup \
        --threads ${small_cpus} \
        --verbosity ${params.bcftools_verbosity} \
        -r ${chr} \
        -f ${params.bcftools_reference} \
        -Ou \
        ${bam} | bcftools call \
        --threads ${small_cpus} \
        --ploidy ${params.bcftools_ploidy} \
        --verbosity ${params.bcftools_verbosity} \
        -mv -Oz -o ${out_vcf}
    """
}


process CONCAT_VCFS {
    
    tag "${sample_id}_Concat_VCFs"

    publishDir "${params.outdir}/03.variant_calling/${sample_id}", mode: 'link', pattern: "*.raw.vcf.gz"

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path("*.raw.vcf.gz"), emit: raw_vcf

    script:
    def out_vcf = "${sample_id}.raw.vcf.gz"
    """
    bcftools concat \
        --threads ${task.cpus} \
        -n \
        -Oz -o ${out_vcf} \
        ${vcfs.join(' ')}
    """
}