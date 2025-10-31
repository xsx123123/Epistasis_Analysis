process BCFTOOLS_MERGE {
    publishDir "${params.outdir}/04_variant_calling/", mode: 'copy', pattern: "merge.sort.vcf.gz*"

    input:
    path vcfs // Channel of all sorted and indexed VCFs

    output:
    path "merge.sort.vcf.gz", emit: merged_vcf
    path "merge.sort.vcf.gz.csi", emit: merged_vcf_csi
    path "merge.sort.vcf.gz.tbi", emit: merged_vcf_tbi

    script:
    """
    ${params.bcftools_path} merge --threads ${task.cpus} \
        ${vcfs.join(' ')} -O z -o merge.vcf.gz
    ${params.bcftools_path} sort \
        merge.vcf.gz -O z -o merge.sort.vcf.gz
    ${params.bcftools_path} index --threads ${task.cpus} \
        -t merge.sort.vcf.gz -o merge.sort.vcf.tbi
    ${params.bcftools_path} index --threads ${task.cpus} \
        -c merge.sort.vcf.gz -o merge.sort.vcf.csi
    """
}