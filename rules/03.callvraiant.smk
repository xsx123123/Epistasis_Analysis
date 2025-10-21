#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
from loguru import logger
# ----- rule ----- #
rule call_variant_by_bcftools:
    input:
        bam = '../02.mapping/bwa_mem2/{sample}.dup.bam',
    output:
        vcf = '../03.call_variant/{sample}.raw.vcf.gz',
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

rule sort_index_bcftools:
    input:
        vcf = '../03.call_variant/{sample}.raw.vcf.gz',
    output:
        sort_vcf = '../03.call_variant/{sample}.sort.vcf.gz',
        sort_vcf_index = '../03.call_variant/{sample}.sort.vcf.csi',
        sort_vcf_tbi_index = '../03.call_variant/{sample}.sort.vcf.tbi',
    message:
        "Running bcftools sort & index on vcf : {input.vcf}",
    log:
        "../logs/03.call_variant/{sample}_variant_sort_index.log",
    benchmark:
            "../benchmarks/{sample}_variant_sort_index_benchmark.txt",
    threads: 
        config["threads"]["bcftools"],
    params:
        path = config["bcftools"]['path'],
    shell:
        """
        ( {params.path} sort --threads {threads} \
          {input.vcf} -O z -o {output.sort_vcf}  && \
          {params.path} index --threads {threads} \
          -t {output.sort_vcf} -o {output.sort_vcf_tbi_index} && \
          {params.path} index --threads {threads} \
          -c {output.sort_vcf} -o {output.sort_vcf_index} ) &>{log}
        """

rule merg_vcf:
    input: 
        vcf_list = expand('../03.call_variant/{sample}.sort.vcf.gz', sample=samples.keys()),
    output:
        merge_vcf = '../03.call_variant/merge.vcf.gz',
        merge_sort_vcf = '../03.call_variant/merge.sort.vcf.gz',
        merge_index_vcf = '../03.call_variant/merge.sort.vcf.gz.csi',
        merge_tbi_index_vcf = '../03.call_variant/merge.sort.vcf.gz.tbi',
    message:
        "Running bcftools merge & sort & index on vcf : {input.vcf}",
    log:
        "../logs/03.call_variant/merge_variant_sort_index.log",
    benchmark:
            "../benchmarks/merge_variant_sort_index_benchmark.txt",
    threads: 
        config["threads"]["bcftools"],
    params:
        path = config["bcftools"]['path'],
    shell:
        """
        ( {params.path} merge --threads {threads} \
          {input.vcf_list} -O z -o {output.merge_vcf}  && \
          {params.path} sort --threads {threads} \
          {output.merge_vcf} -O z -o {output.merge_sort_vcf} && \
          {params.path} index --threads {threads} \
          -t {output.merge_sort_vcf} -o {output.merge_tbi_index_vcf} && \
          {params.path} index --threads {threads} \
          -c {output.merge_sort_vcf} -o {output.merge_index_vcf} ) &>{log}
        """

rule merg_vcf_filter:
    input: 
        merge_vcf = '../03.call_variant/merge.sort.vcf.gz',
    output:
        filter_vcf = '../03.call_variant/merge_filter.sort.vcf.gz',
        filter_index_vcf = '../03.call_variant/merge_filter.sort.vcf.gz.csi',
        filter_tbi_index_vcf = '../03.call_variant/merge_filter.sort.vcf.gz.tbi',
    message:
        "Running bcftools filter  F_MISSING <= 0.1 && MAF >= 0.05 : {input.vcf}",
    log:
        "../logs/03.call_variant/merge_variant_filter.log",
    benchmark:
            "../benchmarks/merge_variant_filter_benchmark.txt",
    threads: 
        config["threads"]["bcftools"],
    params:
        path = config["bcftools"]['path'],
    shell:
        """
        ( {params.path} view --threads {threads} \
          -v snps,indels \
          -i 'F_MISSING <= 0.1 && MAF >= 0.05' \
          {input.merge_vcf} -O z -o {output.filter_vcf} && \
          {params.path} index --threads {threads} \
          -t {output.filter_vcf} -o {output.filter_tbi_index_vcf} && \
          {params.path} index --threads {threads} \
          -c {output.filter_vcf} -o {output.filter_index_vcf}
          ) &>{log}
        """

rule SnpEff_annotation:
    input:
        filter_vcf = '../03.call_variant/merge_filter.sort.vcf.gz',
    output:
        annotated_csv = '../03.call_variant/merge_filter.sort.annotation.csv',
        annotated_html = '../03.call_variant/merge_filter.sort.annotation.html',
        annotated_vcf = '../03.call_variant/merge_filter.sort.annotation.vcf',
    message:
        "Running filter vcf annotation : {input.filter_vcf}",
    log:
        "../logs/03.call_variant/merge_variant_filter_annotation.log",
    benchmark:
            "../benchmarks/merge_variant_filter_annotation_benchmark.txt",
    conda:
        "../envs/snpEff.yaml",
    threads: 
        config["threads"]["snpEff"],
    params:
        config = config["snpEff"]['config'],
        genome_name = config["snpEff"]['genome_name'],
    shell:
        """
        snpEff  -csvStats {output.annotated_csv} \
                -s {output.annotated_html}  \
                -c {params.config} -v \
                -ud 500 {params.genome_name} \
                {input.filter_vcf} > {output.annotated_vcf} 2>{log}
        """

rule variant_stats:
    input:
        vcf = '../03.call_variant/{sample}.sort.vcf.gz',
    output:
        stats = '../03.call_variant/variant_stats/{sample}.vcf.stats.txt',
    message:
        "Running bcftools stats on call variant : {input.vcf}",
    log:
        "../logs/03.call_variant/{sample}_variant_stats.log",
    benchmark:
            "../benchmarks/{sample}_variant_stats_benchmark.txt",
    threads: 
        config["threads"]["bcftools"],
    params:
        reference = config["bcftools"]['reference'],
        path = config["bcftools"]['path'],
    shell:
        """
        {params.path} stats --threads {threads} \
                            -f {params.reference} \
                            {input.vcf}  > {output.stats} 2>{log}
        """

rule variant_stats_multiqc:
    input:
        fastqc_files_r1 = expand("../03.call_variant/variant_stats/{sample}.vcf.stats.txt", sample=samples.keys()),
    output:
        report_dir = directory("../03.call_variant/variant_stats_multiqc/")
    conda:
        "../envs/multiqc.yaml",
    message:
        "Running MultiQC to bcftools variant stats",
    params:
        bcftools_reports = "../03.call_variant/variant_stats/",
        report = "multiqc_bcftools_report.html",
        title = "bcftools-multiqc-report",
    log:
        "../logs/03.call_variant/multiqc-variant-stats.log",
    benchmark:
        "../benchmarks/multiqc-variant-stats_benchmark.txt",
    shell:
        """
        multiqc {params.bcftools_reports} \
                --outdir {output.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """

# ----- rule ----- #