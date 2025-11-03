// 02.qualimap_qc.nf
process QUALIMAP_QC {
    tag "${sample_id}_QUALIMAP_QC"
    
    publishDir "${params.outdir}/02.mapping/${sample_id}", mode: 'link'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}_qualimap/qualimapReport.html", emit: html_report
    path "${sample_id}_qualimap/genome_results.txt", emit: txt_report
    path "${sample_id}_qualimap/raw_data_qualimapReport", emit: raw_data  // 用于 MultiQC

    script:
    """
    qualimap bamqc \
             -nt ${task.cpus} \
             -bam ${bam} \
             -gff ${params.qualimap_genome_gff} \
             -outdir ${sample_id}_qualimap \
             -outformat ${params.qualimap_format} \
             --java-mem-size=${params.qualimap_mem}
    """
}