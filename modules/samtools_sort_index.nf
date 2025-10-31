process SAMTOOLS_SORT_INDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/03_mapping/${sample_id}", mode: 'copy', pattern: "*.sort.bam*"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("*.sort.bam"), path("*.sort.bam.bai"), emit: sorted_indexed_bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sort.bam ${bam}
    samtools index -@ ${task.cpus} ${sample_id}.sort.bam
    """
}
