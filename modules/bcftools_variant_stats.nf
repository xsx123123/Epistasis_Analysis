process BCFTOOLS_VARIANT_STATS {
    tag "$sample_id"
    publishDir "${params.outdir}/04_variant_calling/variant_stats/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(csi), path(tbi)

    output:
    path "*.vcf.stats.txt", emit: variant_stats_report

    conda "${baseDir}/envs/bwa2.yaml" // bcftools is in bwa2.yaml

    script:
    """
    ${params.bcftools_path} stats \
        --threads ${task.cpus} \
        -f ${params.bcftools_reference} \
        ${vcf} > ${sample_id}.vcf.stats.txt
    """
}