process MULTIQC {
    publishDir "${params.outdir}/${output_dir}", mode: 'copy'

    input:
    path report_filesz
    val output_dir
    val report_name
    val report_title

    output:
    path "*.html", emit: multiqc_report

    conda "${baseDir}/envs/multiqc.yaml"

    script:
    """
    multiqc ${report_files.join(' ')} \
            --outdir . \
            -i "${report_title}" \
            -n "${report_name}"
    """
}