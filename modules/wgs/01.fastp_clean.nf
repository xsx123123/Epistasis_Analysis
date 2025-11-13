// 01.fastp_clean.nf
process FASTP_CLEAN {
    tag "${sample_id}_FASTP_CLEAN"
    
    publishDir "${params.outdir}/01.fastp_clean", mode: 'link'

    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}_clean_R1.fq.gz"), path("${sample_id}_clean_R2.fq.gz"), emit: trimmed_reads
    path("${sample_id}.fastp.json"), emit: json
    path("${sample_id}.html")      , emit: html

    script:
    """
    fastp -i ${read1} -I ${read2} \
          -o ${sample_id}_clean_R1.fq.gz -O ${sample_id}_clean_R2.fq.gz \
          -j ${sample_id}.fastp.json -h ${sample_id}.html \
          --thread ${task.cpus}
    """
}