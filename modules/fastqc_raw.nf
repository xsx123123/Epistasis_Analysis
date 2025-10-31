process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1, stageAs: 'read1.fq.gz'), path(read2, stageAs: 'read2.fq.gz')
    
    output:
    path("*.html")
    path("*.zip")
    
    script:
    """
    fastqc -t ${task.cpus} ${params.raw_data_path}/${sample_id}_R1.fq.gz ${params.raw_data_path}/${sample_id}_R2.fq.gz
    """
}
