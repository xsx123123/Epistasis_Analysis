// 03.multiqc_variant_stats.nf
process MULTIQC_VARIANT_STATS {
    tag 'Variant_STATS'
    
    publishDir "${params.outdir}/03.variant_calling/variant_stats_multiqc/", mode: 'link'

    input:
    path stats_files

    output:
    path "variant_multiqc_report.html", emit: multiqc_report

    script:
    """
    multiqc ${stats_files.join(' ')} \
            --outdir . \
            -i "BCFtools Variant Stats MultiQC Report" \
            -n variant_multiqc_report.html
    """
}
