// 02.bwa_index.nf 
process BWA_INDEX {
    tag "${genome}_bwa-mem2_index"

    publishDir "${params.outdir}/02.reference", mode: 'link'

    input:
    path genome

    output:
    path "Lsat_Salinas_v11.genome.fasta*", emit: index 

    script:
    """
    ${params.bwa_mem2_path} index ${genome}
    """
}
