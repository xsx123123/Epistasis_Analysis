#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ----- rule ----- #
rule check_md5:
    input:
        md5 = expand(os.path.join(config['raw_data_path'], '{sample}'),
             sample=samples.keys()),
    output:
        md5_check = "../01.qc/md5_check.tsv",
    message:
        "Running FastQC on {input.md5}",
    params:
        checkmd5_name = config['checkmd5']['name'],
        log_file = "../logs/01.qc/md5_check.log",
    threads: 
        config['threads']['md5_check']
    shell:
        """
        ./scripts/md5_checker_rs/target/release/md5_checker_rs {input.md5} \
                -f {params.checkmd5_name} \
                -t {threads} \
                -o {output.md5_check} \
                --log-file {params.log_file}
        """

rule short_read_qc_r1:
    input:
        md5_check = "../01.qc/md5_check.tsv",
        r1 = os.path.join(config["raw_data_path"],"{sample}", "{sample}" + config['r1_suffix']),
    output:
        r1_html = "../01.qc/short_read_qc_r1/{sample}."+config['r1_suffix']+".html",
        r1_zip = "../01.qc/short_read_qc_r1/{sample}."+config['r1_suffix']+".zip",
    conda:
        "../envs/fastqc.yaml",
    log:
        r1 = "../logs/01.short_read_qc_r1/{sample}.r1.fastqc.log",
    params:
        out_dir = "../01.qc/short_read_qc_r1/",
    message:
        "Running FastQC on {input.r1}",
    threads: 1
    shell:
        """
        fastqc {input.r1} \
               -o {params.out_dir} \
               --threads {threads} &> {log.r1}
        """

rule short_read_qc_r2:
    input:
        md5_check = "../01.qc/md5_check.tsv",
        r2 = os.path.join(config["raw_data_path"],"{sample}", "{sample}" + config['r2_suffix']),
    output:
        r2_html = "../01.qc/short_read_qc_r1/{sample}."+config['r2_suffix']+".html",
        r2_zip = "../01.qc/short_read_qc_r1/{sample}."+config['r2_suffix']+".zip",
    conda:
        "../envs/fastqc.yaml",
    log:
        r2 = "../logs/01.short_read_qc_r2/{sample}.r2.fastqc.log",
    params:
        out_dir = "../01.qc/short_read_qc_r2",
    message:
        "Running FastQC on {input.r2}",
    threads: 1
    shell:
        """
        fastqc {input.r2} \
               -o {params.out_dir} \
               --threads {threads} &> {log.r2}
        """

# logger.info('Run MultiQC to summarize R1 fastqc QC reports')
rule short_read_multiqc_r1:
    input:
        fastqc_files_r1 = expand("../01.qc/short_read_qc_r1/{sample}."+config['r1_suffix']+".zip", sample=samples.keys()),
    output:
        report_dir = directory("../01.qc/short_read_r1_multiqc/")
    conda:
        "../envs/multiqc.yaml",
    message:
        "Running MultiQC to aggregate R1 FastQC reports",
    params:
        fastqc_reports = "../01.qc/short_read_qc_r1",
        report = "multiqc_r1_raw-data_report.html",
        title = "r1-raw-data-multiqc-report",
    log:
        "../logs/01.multiqc/multiqc-r1.log",
    shell:
        """
        multiqc {params.fastqc_reports} \
                --outdir {output.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """

# logger.info('Run MultiQC to summarize R2 fastqc QC reports')
rule short_read_multiqc_r2:
    input:
        fastqc_files_r2 = expand("../01.qc/short_read_qc_r1/{sample}."+config['r2_suffix']+".zip", sample=samples.keys()),
    output:
        report_dir = directory("../01.qc/short_read_r2_multiqc/")
    conda:
        "../envs/multiqc.yaml",
    message:
        "Running MultiQC to aggregate R2 FastQC reports",
    params:
        fastqc_reports = "../01.qc/short_read_qc_r2",
        report = "multiqc_r2_raw-data_report.html",
        title = "r2-raw-data-multiqc-report",
    log:
        "../logs/01.multiqc/multiqc-r2.log",
    shell:
        """
        multiqc {params.fastqc_reports} \
                --outdir {output.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """
# ----- rule ----- #