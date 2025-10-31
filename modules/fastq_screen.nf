process FASTQ_SCREEN {
    tag "$sample_id"
    publishDir "${params.outdir}/01_fastq_screen/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read1, stageAs: '${sample_id}_R1.fq.gz'), path(read2, stageAs: '${sample_id}_R2.fq.gz')

    output:
    path "*_screen.txt", emit: screen_txt
    path "*_screen.html", emit: screen_html

    conda "${baseDir}/envs/fastqc.yaml"

    script:
    """
    ${params.fastq_screen_path} --threads ${task.cpus} \
                                --force \
                                --subset ${params.fastq_screen_subset} \
                                --aligner ${params.fastq_screen_aligner} \
                                --conf ${params.fastq_screen_conf} \
                                --outdir . \
                                ${params.raw_data_path}/${sample_id}_R1.fq.gz
    ${params.fastq_screen_path} --threads ${task.cpus} \
                                --force \
                                --subset ${params.fastq_screen_subset} \
                                --aligner ${params.fastq_screen_aligner} \
                                --conf ${params.fastq_screen_conf} \
                                --outdir . \
                                ${params.raw_data_path}/${sample_id}_R2.fq.gz
    """
}