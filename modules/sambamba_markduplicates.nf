process SAMBAMBA_MARKDUPLICATES {
    tag "$sample_id"
    publishDir "${params.outdir}/03_mapping/${sample_id}", mode: 'copy', pattern: "*.dup.bam"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*.dup.bam"), emit: marked_duplicates_bam

    script:
    """
    ${params.sambamba_path} markdup \
        --nthreads ${task.cpus} \
        --show-progress \
        ${bam} \
        ${sample_id}.dup.bam
    """
}
