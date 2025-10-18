#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
from loguru import logger
# ----- rule ----- #
if judge_bwa_index(config = config) == True:
    logger.info(f"Epistasis Analysis pipeline will rebuild BWA index")
    rule build_bwa_index:
        input:
            genome = config["genome"],
        output:
            config["bwa_mem2"]["index"] + '.0123',
            config["bwa_mem2"]["index"] + '.amb',
            config["bwa_mem2"]["index"] + '.ann',
            config["bwa_mem2"]["index"] + '.bwt.2bit.64',
            config["bwa_mem2"]["index"] + '.pac',
            config["bwa_mem2"]["index"] + '.alt',
        conda:
            "../envs/bwa2.yaml", 
        log:
            "logs/02.mapping/bwa_mem2.log"
        message:
            "Building BWA-mem2 index for {input.genome}",
        benchmark:
            "../benchmarks/BWA-mem2_index_benchmark.txt",
        shell:
            """
            bwa-mem2 index \
                     -p {config["bwa_mem2"]["index"]} \
                      {input.genome} 2>{log}
            """

rule bwa_mapping:
    input:
        r1 = "../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz",
        r2 = "../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz",
    output:
        bam = '../02.mapping/bwa_mem2/{sample}.bam',
    conda:
        "../envs/bwa2.yaml",
    log:
        "../logs/02.mapping/bwa_mem2_{sample}.log",
    message:
        "Running bwa-mem2 mapping on {input.r1} and {input.r2}",
    benchmark:
            "../benchmarks/{sample}_bwa_mem2_benchmark.txt",
    params:
        index = config["bwa_mem2"]["index"],
        pl = config["bwa_mem2"]["PL"],
    threads: 
        config["threads"]["bwa_mem2"],
    shell:
        """
        ( bwa-mem2 mem -t {threads} \
        -R '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:{params.pl}\tLB:{wildcards.sample}' \
        {params.index} \
        {input.r1} {input.r2} | samtools view -@ {threads} \
        -Sbh -o {output.bam} ) 2>{log}
        """

rule sort_index:
    input:
        bam = '../02.mapping/bwa_mem2/{sample}.bam',
    output:
        sort_bam = '../02.mapping/bwa_mem2/{sample}.sort.bam',
        sort_bam_bai = '../02.mapping/bwa_mem2/{sample}.sort.bam.bai',
    conda:
        "../envs/bwa2.yaml",
    message:
        "Running samtools sort & index for {input.bam}",
    log:
        "../logs/02.mapping/bwa_sort_index_{sample}.log",
    benchmark:
            "../benchmarks/{sample}_bam_sort_index_benchmark.txt",
    threads: 
        config["threads"]["samtools"],
    shell:
        """
        (samtools sort -@ {threads} -o {output.sort_bam} {input.bam} &&
        samtools index -@ {threads} {output.sort_bam}) 2>{log}
        """

rule sambamba_MarkDuplicates:
    input:
        bam = '../02.mapping/bwa_mem2/{sample}.sort.bam',
    output:
        duplicates_bam = '../02.mapping/bwa_mem2/{sample}.dup.bam',
    log:
        "../logs/02.mapping/duplicates_{sample}.log",
    message:
        "Running BAM MarkDuplicates for {input.bam}",
    benchmark:
            "../benchmarks/{sample}_bam_MarkDuplicates_benchmark.txt",
    params:
        sambamba = config['software']['sambamba']
    threads: 8
    shell:
        """
        {params.sambamba} markdup \
                          --nthreads=NTHREADS {threads} \
                          --show-progress \
                          {input.bam} \
                          {output.duplicates_bam} 2>{log}
        """

rule Duplicates_bam_index:
    input:
        bam = '../02.mapping/bwa_mem2/{sample}.dup.bam',
    output:
        duplicates_bam_bai = '../02.mapping/bwa_mem2/{sample}.dup.bam.bai',
    conda:
        "../envs/bwa2.yaml",
    message:
        "Running build index for MarkDuplicates of BAM {input.bam}",
    log:
        "../logs/02.mapping/Duplicates_bam_index_{sample}.log",
    benchmark:
            "../benchmarks/{sample}_Dup_bam_index_benchmark.txt",
    threads: 
        config["threads"]["samtools"],
    shell:
        """
        samtools index -@ {threads} {input.bam} 2>{log}
        """

rule bam_coverage:
    input:
        bam = '../02.mapping/bwa_mem2/{sample}.dup.bam',
    output:
        dist = "../02.mapping/mosdepth_coverage/{sample}.mosdepth.global.dist.txt",
        summary = "../02.mapping/mosdepth_coverage/{sample}.mosdepth.summary.txt",
    conda:
        "../envs/mosdepth.yaml",
    message:
        "Running Coverage for MarkDuplicates of BAM : {input.bam}",
    log:
        "../logs/02.mapping/bam_coverage_{sample}.log",
    benchmark:
            "../benchmarks/{sample}_Dup_bam_coverage_benchmark.txt",
    params:
        prefix = "../02.mapping/mosdepth_coverage/{sample}"
    threads: 
        config["threads"]["qualimap"],
    shell:
        """
        mosdepth -n \
                 --fast-mode \
                 -t {threads} \
                 --by 500 \
                 {params.prefix} \
                 {input.bam} 2>{log}
        """

rule qualimap_qc:
    input:
        bam = '../02.mapping/bwa_mem2/{sample}.dup.bam',
    output:
        qualimap_report_html = '../02.mapping/qualimap_report/{sample}_qualimap_report.html',
    conda:
        "../envs/qualimap.yaml",
    message:
        "Running qualimap qc for MarkDuplicates of BAM : {input.bam}",
    log:
        "../logs/02.mapping/qualimap_report_{sample}.log",
    benchmark:
            "../benchmarks/{sample}_Dup_bam_qualimap_benchmark.txt",
    params:
        genome_gff = config['qualimap']["genome_gff"],
        outformat = config['qualimap']["format"],
        prefix_dir = '../02.mapping/qualimap_report/',
    threads: 
        config["threads"]["qualimap"],
    shell:
        """
        qualimap bamqc \
                 -nt {threads} \
                 -bam {input.bam} \
                 -gff {params.genome_gff} \
                 -outdir {params.prefix_dir} \
                 -outfile {output.qualimap_report_html} \
                 -outformat {params.outformat} 2>{log}
        """

rule samtools_flagst:
    input:
        bam = '../02.mapping/bwa_mem2/{sample}.dup.bam',
    output:
        samtools_flagstat = '../02.mapping/samtools_flagstat/{sample}_dup_bam_flagstat.tsv',
    conda:
        "../envs/bwa2.yaml",
    message:
        "Running flagst for MarkDuplicates of BAM : {input.bam}",
    log:
        "../logs/02.mapping/bam_dup_lagstat_{sample}.log",
    benchmark:
            "../benchmarks/{sample}_Dup_bam_lagstat_benchmark.txt",
    threads: 
        config["threads"]["samtools_flagstat"],
    shell:
        """
        samtools flagstat \
                 -@ {threads} \
                 -O tsv \
                 {input.bam} > {output.samtools_flagstat} 2>{log}
        """

rule samtools_stats:
    input:
        bam = '../02.mapping/bwa_mem2/{sample}.dup.bam',
    output:
        samtools_stats = '../02.mapping/samtools_stats/{sample}_dup_bam_stats.tsv',
    conda:
        "../envs/bwa2.yaml",
    message:
        "Running stats for MarkDuplicates of BAM : {input.bam}",
    log:
        "../logs/02.mapping/bam_dup_stats_{sample}.log",
    benchmark:
            "../benchmarks/{sample}_Dup_bam_stats_benchmark.txt",
    threads: 
        config["threads"]["samtools_stats"],
    params:
        reference = config["samtools_stats"]['reference']
    shell:
        """
        samtools stats \
                 -@ {threads} \
                 --reference {params.reference} \
                 {input.bam} > {output.samtools_stats} 2>{log}
        """
# ----- rule ----- #