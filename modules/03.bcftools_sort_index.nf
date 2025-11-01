// 03.bcftools_sort_index.nf 
process BCFTOOLS_SORT_INDEX {
    tag "${sample_id}_BCFTOOLS_SORT_INDEX"

    publishDir "${params.outdir}/03.variant_calling/${sample_id}", mode: 'link', pattern: "*.sort.vcf.gz*"
    
    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), 
          path("${sample_id}.sort.vcf.gz"), 
          path("${sample_id}.sort.vcf.gz.csi"), 
          path("${sample_id}.sort.vcf.gz.tbi"), 
    emit: sorted_indexed_vcf

    script:
    """
    bcftools sort \
        ${vcf} -O z -o ${sample_id}.sort.vcf.gz
    
    bcftools index -f -t \
        --threads ${task.cpus} \
        ${sample_id}.sort.vcf.gz
    
    bcftools index -f \
        --threads ${task.cpus} \
        -c ${sample_id}.sort.vcf.gz
    """
}
