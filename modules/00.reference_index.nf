process SAMTOOLS_FAIDX {
    
    tag "samtools faidx ${fasta.baseName}"

    publishDir "${params.outdir}/reference_index", mode: 'link', pattern: "*.fai"

    input:
    path(fasta)

    output:
    path("${fasta.name}.fai"), emit: fai

    script:
    """
    ln -s ${fasta} .
    samtools faidx ${fasta.name}
    """
}