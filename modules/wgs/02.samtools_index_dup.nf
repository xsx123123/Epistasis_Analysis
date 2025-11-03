// 02.samtools_index_dup.nf
process SAMTOOLS_INDEX_DUP {
    tag "${sample_id}_INDEX_DUP_BAM"

    publishDir "${params.outdir}/02.mapping/${sample_id}", mode: 'link', pattern: "*.dup.bam.bai"
    
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id),path("*.dup.bam.bai"), emit: indexed_dup_bam

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}
