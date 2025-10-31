process FASTP_CLEAN {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp_clean", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}_clean_R1.fq.gz"), path("${sample_id}_clean_R2.fq.gz"), emit: trimmed_reads
    path("*.fastp.json") 
    path("*.html") 
    
    script:
    """
    fastp -i ${params.raw_data_path}/$read1 -I ${params.raw_data_path}/$read2 \
          -o ${sample_id}_clean_R1.fq.gz -O ${sample_id}_clean_R2.fq.gz \
          -j ${sample_id}.fastp.json -h ${sample_id}.html \
          --thread ${task.cpus}
    """
}