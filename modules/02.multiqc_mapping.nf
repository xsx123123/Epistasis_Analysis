// 02.multiqc_mapping.nf
process MULTIQC_MAPPING {
    tag 'Raw_data_MULTIQC'
    publishDir "${params.outdir}/02.mapping/", mode: 'link'

    input:
    path fastqc_files

    output:
    path "mapping_report.html"

    script:
    """
    multiqc . --force --filename mapping_report.html
    """
}