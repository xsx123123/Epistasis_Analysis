process BWA_MEM {
    tag "$sample_id"
    publishDir "${params.outdir}/03_mapping/${sample_id}", mode: 'copy', pattern: "*.bam"

    input:
    tuple val(sample_id), path(fq1), path(fq2)
    path bwa_index_dir

    output:
    tuple val(sample_id), path("*.bam"), emit: bam

    script:
    def read_group = "'@RG\tID:${sample_id}\tSM:${sample_id}\tPL:${params.bwa_mem2_pl}\tLB:${sample_id}'"
    """
    bwa-mem2 mem -t ${task.cpus} \
        -R ${read_group} \
        ${bwa_index_dir}/genome \
        ${fq1} ${fq2} | \
    samtools view -@ ${task.cpus} -Sbh -o ${sample_id}.bam
    """
}
