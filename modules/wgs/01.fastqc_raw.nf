// 01.fastqc_raw.nf
process FASTQC_RAW {
    tag "${sample_id}_RAW_FASTQC"
    
    publishDir "${params.outdir}/01.fastqc_raw/${sample_id}", mode: 'link', pattern: "${sample_id}*.{html,zip}" 

    input:
    tuple val(sample_id), path(read1_file), path(read2_file)
    
    output:
    path("${sample_id}_R1_fastqc.html"), emit: html_r1
    path("${sample_id}_R2_fastqc.html"), emit: html_r2
    path("${sample_id}_R1_fastqc.zip"), emit: zip_r1
    path("${sample_id}_R2_fastqc.zip"), emit: zip_r2
    
    script:
    """
    fastqc -t ${task.cpus} ${read1_file} ${read2_file} -o . 
    """
}