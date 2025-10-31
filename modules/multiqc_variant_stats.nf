process MULTIQC_VARIANT_STATS {
    publishDir "${params.outdir}/04_variant_calling/variant_stats_multiqc/", mode: 'copy'

    input:
    path stats_files // Channel of all variant stats files

    output:
    path "multiqc_report.html", emit: multiqc_report

    conda "${baseDir}/envs/multiqc.yaml"

    script:
    """
    multiqc ${stats_files.join(' ')} \
            --outdir . \
            -i "BCFtools Variant Stats MultiQC Report" \
            -n multiqc_report.html
    """
}
