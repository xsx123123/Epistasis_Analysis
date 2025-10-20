#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ------- rule ------- #
rule seq_preprocessor:
    input:
        md5 = expand(os.path.join(config['raw_data_path'], '{sample}'),
                     sample=samples.keys()),
    output:
        md5_check = os.path.join(config['raw_data_path'],config['convert_md5']),
        md5_check_json = os.path.join(config['raw_data_path'],config['convert_md5'],"raw_data_md5.json"),
    message:
        "Running seq_preprocessor on raw data data",
    benchmark:
        "../benchmarks/seq_preprocessor.txt",
    params:
        raw_data_path = config['raw_data_path'],
        md5 = config['md5']
    threads: 1
    shell:
        """
        ./scripts/seq_preprocessor/target/release/seq_preprocessor -i  {params.raw_data} \
                -o {output.md5_check} \
                --md5-name {params.md5} \              
                --json-report {output.md5_check_json} &>{log}
        """

rule check_md5:
    input:
        md5_check_json = os.path.join(config['raw_data_path'],config['convert_md5'],"raw_data_md5.json"),
    output:
        md5_check = "../01.qc/md5_check.tsv",
    message:
        "Running md5 check on raw data files on {input.md5} directory",
    benchmark:
        "../benchmarks/md5_check_benchmark.txt",
    params:
        md5_check = os.path.join(config['raw_data_path'],config['convert_md5']),
        log_file = "../logs/01.qc/md5_check.log",
    threads: 
        config['threads']['md5_check']
    shell:
        """
        ./scripts/json_md5_verifier/target/release/json_md5_verifier -t  {threads} \
                -i {input.md5_check_json} \
                -b {params.md5_check} \              
                -o {output.md5_check} \
                --log-file {params.log_file}
        """