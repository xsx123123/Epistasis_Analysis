// 03.bcftools_filter.nf 
process BCFTOOLS_FILTER {
    tag "Merge_Bcftools_filter_variant"

    publishDir "${params.outdir}/03.variant_calling/", mode: 'link', pattern: "merge_filter.sort.vcf.gz*"
    
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
    bcftools view --threads ${task.cpus} \
        -v snps,indels \
        -i 'F_MISSING <= 0.1 && MAF >= 0.05' \
        ${merged_vcf} -O z -o merge_filter.sort.vcf.gz
    bcftools index --threads ${task.cpus} \
        -t merge_filter.sort.vcf.gz -o merge_filter.sort.vcf.gz.tbi
    bcftools index --threads ${task.cpus} \
        -c merge_filter.sort.vcf.gz -o merge_filter.sort.vcf.gz.csi
    """
}