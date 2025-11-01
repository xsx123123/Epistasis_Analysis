// 01.fastq_screen.nf
process FASTQ_SCREEN {
    tag "${sample_id}_FASTQ_SCREEN"
    
    publishDir "${params.outdir}/01.fastq_screen/", mode: 'link'
        
    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "*_screen.txt", emit: screen_txt
    path "*_screen.html", emit: screen_html

    script:
    """
    fastq_screen --threads ${task.cpus} \
                                --force \
                                --subset ${params.fastq_screen_subset} \
                                --aligner ${params.fastq_screen_aligner} \
                                --conf ${params.fastq_screen_conf} \
                                --outdir . \
                                ${params.raw_data_path}/${sample_id}_R1.fq.gz
    fastq_screen --threads ${task.cpus} \
                                --force \
                                --subset ${params.fastq_screen_subset} \
                                --aligner ${params.fastq_screen_aligner} \
                                --conf ${params.fastq_screen_conf} \
                                --outdir . \
                                ${params.raw_data_path}/${sample_id}_R2.fq.gz
    """
}