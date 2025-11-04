// 02.bam_coverage.nf 
process BAM_COVERAGE {
    
    tag "${sample_id}_bam_coverage"

    publishDir "${params.outdir}/02.mapping/${sample_id}", mode: 'link'
    
    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "*.mosdepth.global.dist.txt", emit: global_dist
    path "*.mosdepth.summary.txt", emit: summary

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