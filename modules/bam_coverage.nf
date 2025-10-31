process BAM_COVERAGE {
    tag "$sample_id"
    publishDir "${params.outdir}/03_mapping/mosdepth_coverage/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "*.mosdepth.global.dist.txt", emit: global_dist
    path "*.mosdepth.summary.txt", emit: summary

    conda "${baseDir}/envs/mosdepth.yaml"

    script:
    """
    mosdepth -n \
             --fast-mode \
             -t ${task.cpus} \
             --by 500 \
             ${sample_id} \
             ${bam}
    """
}