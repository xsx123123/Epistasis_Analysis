// 01.fastqc_raw.nf 
process FASTQC_RAW {
    tag "${sample_id}_RAW_FASTQC"
    
    publishDir "${params.outdir}/01.fastqc_raw",mode: 'link',pattern: '*.{html,zip}'

    input:
    tuple val(sample_id),path(read1_file),path(read2_file)
    
    output:
    path("*.html"), emit: html
    path("*.zip"), emit: zip
    
    script:
    """
    fastqc -t ${task.cpus} ${read1_file} -o . && \
    fastqc -t ${task.cpus} ${read2_file} -o .
    """
}
