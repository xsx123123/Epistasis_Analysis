process SAMTOOLS_INDEX_DUP {
    tag "$sample_id"
    publishDir "${params.outdir}/03_mapping/${sample_id}", mode: 'copy', pattern: "*.dup.bam.bai"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("*.dup.bam"), path("*.dup.bam.bai"), emit: indexed_dup_bam

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}
