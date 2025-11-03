// 01.multiqc.nf
process MULTIQC {
    tag 'Raw_data_MULTIQC'
    
    publishDir "${params.outdir}/01.multiqc", mode: 'link'

    input:
    path fastqc_files

    output:
    path "raw_clean_data_report.html"

    script:
    """
    multiqc . --force --filename raw_clean_data_report.html
    """
}