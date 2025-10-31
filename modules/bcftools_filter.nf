process BCFTOOLS_FILTER {
    publishDir "${params.outdir}/04_variant_calling/", mode: 'copy', pattern: "merge_filter.sort.vcf.gz*"

    input:
    path merged_vcf
    path merged_vcf_csi
    path merged_vcf_tbi

    output:
    path "merge_filter.sort.vcf.gz", emit: filtered_vcf
    path "merge_filter.sort.vcf.gz.csi", emit: filtered_vcf_csi
    path "merge_filter.sort.vcf.gz.tbi", emit: filtered_vcf_tbi

    script:
    """
    ${params.bcftools_path} view --threads ${task.cpus} \
        -v snps,indels \
        -i 'F_MISSING <= 0.1 && MAF >= 0.05' \
        ${merged_vcf} -O z -o merge_filter.sort.vcf.gz
    ${params.bcftools_path} index --threads ${task.cpus} \
        -t merge_filter.sort.vcf.gz -o merge_filter.sort.vcf.gz.tbi
    ${params.bcftools_path} index --threads ${task.cpus} \
        -c merge_filter.sort.vcf.gz -o merge_filter.sort.vcf.gz.csi
    """
}