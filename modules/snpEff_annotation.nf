process SNPEFF_ANNOTATION {
    publishDir "${params.outdir}/04_variant_calling/", mode: 'copy', pattern: "merge_filter.sort.annotation.*"

    input:
    path filtered_vcf
    path filtered_vcf_csi
    path filtered_vcf_tbi

    output:
    path "merge_filter.sort.annotation.csv", emit: annotated_csv
    path "merge_filter.sort.annotation.html", emit: annotated_html
    path "merge_filter.sort.annotation.vcf", emit: annotated_vcf

    conda "${baseDir}/envs/snpEff.yaml"

    script:
    """
    snpEff -csvStats merge_filter.sort.annotation.csv \
           -s merge_filter.sort.annotation.html \
           -c ${params.snpEff_config} -v \
           -ud 500 ${params.snpEff_genome_name} \
           ${filtered_vcf} > merge_filter.sort.annotation.vcf
    """
}

