process BCFTOOLS_SORT_INDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/04_variant_calling/${sample_id}", mode: 'copy', pattern: "*.sort.vcf.gz*"

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("*.sort.vcf.gz"), path("*.sort.vcf.csi"), path("*.sort.vcf.tbi"), emit: sorted_indexed_vcf

    script:
    """
    ${params.bcftools_path} sort \
        ${vcf} -O z -o ${sample_id}.sort.vcf.gz
    ${params.bcftools_path} index --threads ${task.cpus} \
        -t ${sample_id}.sort.vcf.gz -o ${sample_id}.sort.vcf.tbi
    ${params.bcftools_path} index --threads ${task.cpus} \
        -c ${sample_id}.sort.vcf.gz -o ${sample_id}.sort.vcf.csi
    """
}
