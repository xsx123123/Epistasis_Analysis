// 01.fastq_screen.nf
process FASTQ_SCREEN {
    tag "${sample_id}_FASTQ_SCREEN"
    
    publishDir "${params.outdir}/01.fastq_screen/${sample_id}", mode: 'link'
        
    input:
    tuple val(sample_id), path(read1), path(read2) 

    output:
    path("${sample_id}_R1_screen.txt"), emit: R1_screen_txt
    path("${sample_id}_R1_screen.html"), emit: R1_screen_html
    path("${sample_id}_R2_screen.txt"), emit: R2_screen_txt
    path("${sample_id}_R2_screen.html"), emit: R2_screen_html


    script:
    """
    # Read 1 
    fastq_screen --threads ${task.cpus} \
                 --force \
                 --subset ${params.fastq_screen_subset} \
                 --aligner ${params.fastq_screen_aligner} \
                 --conf ${params.fastq_screen_conf} \
                 --outdir . \
                 ${read1} 
    # Read 2  
    fastq_screen --threads ${task.cpus} \
                 --force \
                 --subset ${params.fastq_screen_subset} \
                 --aligner ${params.fastq_screen_aligner} \
                 --conf ${params.fastq_screen_conf} \
                 --outdir . \
                 ${read2} 
    """
}