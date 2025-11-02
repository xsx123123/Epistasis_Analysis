// 03.bcftools_merge.nf 
process BCFTOOLS_MERGE {

    tag "Bcftools_merge"
    
    publishDir "${params.outdir}/03.variant_calling/", mode: 'link', pattern: "merge.sort.vcf.gz*"
    
    input:
    path vcfs
    path csis
    path tbis
    
    output:
    path "merge.sort.vcf.gz", emit: merged_vcf
    path "merge.sort.vcf.gz.csi", emit: merged_vcf_csi
    path "merge.sort.vcf.gz.tbi", emit: merged_vcf_tbi

    script:
    """
    bcftools merge --threads ${task.cpus} \
        ${vcfs.join(' ')} -O z -o merge.vcf.gz
    bcftools sort \
        merge.vcf.gz -O z -o merge.sort.vcf.gz
    bcftools index --threads ${task.cpus} \
        -t merge.sort.vcf.gz -o merge.sort.vcf.gz.tbi
    bcftools index --threads ${task.cpus} \
        -c merge.sort.vcf.gz -o merge.sort.vcf.gz.csi
    """
}