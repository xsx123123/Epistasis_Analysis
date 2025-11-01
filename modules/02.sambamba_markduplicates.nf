// 02.sambamba_markduplicates.nf
process SAMBAMBA_MARKDUPLICATES {
    tag "${sample_id}_MARKDUPLICATES"
    
    publishDir "${params.outdir}/02.mapping/${sample_id}", mode: 'link'
    
    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.dup.bam"), path("${sample_id}.dup.bam.bai"), emit: marked_duplicates_bam
    
    script:
    """
    sambamba markdup \
        --nthreads ${task.cpus} \
        --show-progress \
        ${bam} \
        ${sample_id}.dup.bam
    samtools index -@ ${task.cpus} ${sample_id}.dup.bam
    """
}
