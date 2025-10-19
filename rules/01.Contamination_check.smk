#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ----- rule ----- #
rule short_read_fastq_screen_r1:
    input:
        r1 = os.path.join(config["raw_data_path"],"{sample}", "{sample}" + config['r1_suffix']),
    output:
        fastq_screen_result = "../01.qc/fastq_screen/{sample}_1_screen.txt",
    log:
        "../logs/01.short_read_qc_r1/{sample}.r1.fastq_screen.log",
    params:
        out_dir = "../01.qc/fastq_screen/",
        conf = config['fastq_screen']['conf'],
        subset = config['fastq_screen']['subset'],
        aligner = config['fastq_screen']['aligner'],
    message:
        "Running fastq_screen on {input.r1}",
    benchmark:
        "../benchmarks/{sample}_r1_fastq_screen_benchmark.txt",
    threads: 
        config['threads']['fastq_screen'],
    shell:
        """
        fastq_screen --threads  {threads} \
                     --force \
                     --subset  {params.subset} \
                     --aligner  {params.aligner} \
                     --conf {params.conf} \
                     --outdir {params.out_dir} \
                     {input.r1} &> {log}
        """

rule short_read_fastq_screen_r2:
    input:
        r2 = os.path.join(config["raw_data_path"],"{sample}", "{sample}" + config['r2_suffix']),
    output:
        fastq_screen_result = "../01.qc/fastq_screen/{sample}_2_screen.txt",
    log:
        "../logs/01.short_read_qc_r2/{sample}.r2.fastq_screen.log",
    params:
        out_dir = "../01.qc/fastq_screen/",
        conf = config['fastq_screen']['conf'],
        subset = config['fastq_screen']['subset'],
        aligner = config['fastq_screen']['aligner'],
    message:
        "Running fastq_screen on {input.r2}",
    benchmark:
        "../benchmarks/{sample}_r2_fastq_screen_benchmark.txt",
    threads: 
        config['threads']['fastq_screen'],
    shell:
        """
        fastq_screen --threads  {threads} \
                     --force \
                     --subset  {params.subset} \
                     --aligner  {params.aligner} \
                     --conf {params.conf} \
                     --outdir {params.out_dir} \
                     {input.r2} &> {log}
        """
# ----- rule ----- #