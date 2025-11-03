// 02.bwa_mem.nf
process BWA_MEM {
    tag "${sample_id}_BWA_MEM"
    
    publishDir "${params.outdir}/02.mapping/${sample_id}", mode: 'link', pattern: "*.bam"
    
    input:
    tuple val(sample_id), path(fq1), path(fq2)
    each path(bwa_index)

    output:
    tuple val(sample_id), path("*.bam"), emit: bam

    script:
    def read_group = "'@RG\tID:${sample_id}\tSM:${sample_id}\tPL:${params.bwa_mem2_pl}\tLB:${sample_id}'"
    """
    bwa-mem2 mem -t ${task.cpus} \
        -R ${read_group} \
        ${params.index_name} \
        ${fq1} ${fq2} | \
    samtools view -@ ${task.cpus} -Sbh -o ${sample_id}.bam
    """
}
