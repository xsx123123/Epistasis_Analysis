process BWA_INDEX {
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path genome

    output:
    path "index"

    script:
    """
    mkdir index
    bwa-mem2 index ${genome}
    """
}
