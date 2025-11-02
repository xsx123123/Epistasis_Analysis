// 02.samtools_sort_index.nf
process SAMTOOLS_SORT_INDEX {
    tag "${sample_id}_SORT_INDEX_BAM"
    
    publishDir "${params.outdir}/02.mapping/${sample_id}", mode: 'link', pattern: "*.sort.bam*"
    
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
