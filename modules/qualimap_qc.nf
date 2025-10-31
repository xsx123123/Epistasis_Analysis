process QUALIMAP_QC {
    tag "$sample_id"
    publishDir "${params.outdir}/03_mapping/qualimap_report/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "qualimapReport.html", emit: html_report
    path "genome_results.txt", emit: txt_report

    conda "${baseDir}/envs/qualimap.yaml"

    script:
    """
    qualimap bamqc \
             -nt ${task.cpus} \
             -bam ${bam} \
             -gff ${params.qualimap_genome_gff} \
             -outdir . \
             -outformat ${params.qualimap_format} \
             --java-mem-size=${params.qualimap_mem}
    """
}
