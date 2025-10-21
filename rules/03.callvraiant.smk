#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
from loguru import logger
# ----- rule ----- #
rule call_variant_by_bcftools:
    input:
        bam = '../02.mapping/bwa_mem2/{sample}.dup.bam',
    output:
        vcf = '../03.call_variant/{sample}.vcf.gz',
    message:
        "Running bcftools mpileup & bcftools call on the sorted and duplicate marked BAM : {input.bam}",
    log:
        "../logs/03.call_variant/{sample}_variant.log",
    benchmark:
            "../benchmarks/{sample}_variant_benchmark.txt",
    threads: 
        config["threads"]["bcftools"],
    params:
        reference = config["bcftools"]['reference'],
        path = config["bcftools"]['path'],
    shell:
        """
        {params.path} mpileup \
                 --threads {threads} \
                 -f {params.reference} \
                 {input.bam} | {params.path} call \
                 --threads {threads} \
                 -mv -Oz -o {output.vcf} &>{log}
        """
# ----- rule ----- #