// 02.bwa_index.nf 
process BWA_INDEX {
    tag "${genome}_BWA_INDEX"

    publishDir "${params.outdir}/02.reference", mode: 'link'

    input:
    path genome

    output:
    path "index", emit: index 

    script:
    """
    mkdir index
    bwa-mem2 index -p index/${params.bwa_index_prefix} ${genome}
    """
}
