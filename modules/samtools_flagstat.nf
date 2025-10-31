process SAMTOOLS_FLAGSTAT {
    tag "$sample_id"
    publishDir "${params.outdir}/03_mapping/samtools_flagstat/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "*.flagstat.tsv", emit: flagstat_report

    conda "${baseDir}/envs/bwa2.yaml" // samtools is in bwa2.yaml

    script:
    """
    samtools flagstat \
             -@ ${task.cpus} \
             -O tsv \
             ${bam} > ${sample_id}.flagstat.tsv
    """
}
