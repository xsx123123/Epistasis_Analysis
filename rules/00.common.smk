#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# loading packages
from loguru import logger
from pathlib import Path
from typing import Dict, Union
from rich import print as rich_print
# Target rule function
def Epistasis(config:dict = None) -> list:
    """
    Epistasis analysis function. This function performs epistasis analysis on the input configuration and returns a list of results.
    """
    # short-read raw-data qc result
    WGS_Epistasis = [
            "../01.qc/md5_check.tsv",
            "../01.qc/short_read_r1_multiqc/",
            "../01.qc/short_read_r2_multiqc/"    
                    
        ]
    if config['fastq_screen']['run']:
        WGS_Epistasis.extend(expand("../01.qc/fastq_screen_r1/{sample}_R1_screen.txt",
                                          sample=samples.keys()))
        WGS_Epistasis.extend(expand("../01.qc/fastq_screen_r2/{sample}_R2_screen.txt",
                                          sample=samples.keys()))
    # short-read trim & clean result
    WGS_Epistasis.extend(expand("../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz",
                                          sample=samples.keys()))
    WGS_Epistasis.extend(expand("../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz",
                                          sample=samples.keys()))
    WGS_Epistasis.extend(expand("../01.qc/short_read_trim/{sample}.fastp.html",
                                          sample=samples.keys()))
    WGS_Epistasis.extend(expand("../01.qc/short_read_trim/{sample}.fastp.json",
                                          sample=samples.keys()))
    WGS_Epistasis.append("../01.qc/multiqc_short_read_trim/")
    WGS_Epistasis.append("../01.qc/fastq_screen_multiqc_r1/")
    WGS_Epistasis.append("../01.qc/fastq_screen_multiqc_r2/")
    # mapping result
    WGS_Epistasis.extend(expand("../02.mapping/bwa_mem2/{sample}.sort.bam",
                                           sample=samples.keys()))
    WGS_Epistasis.extend(expand("../02.mapping/bwa_mem2/{sample}.sort.bam.bai",
                                           sample=samples.keys()))
    WGS_Epistasis.extend(expand('../02.mapping/bwa_mem2/{sample}.dup.sort.bam',
                                           sample=samples.keys()))
    WGS_Epistasis.extend(expand('../02.mapping/bwa_mem2/{sample}.dup.sort.bam.bai',
                                           sample=samples.keys()))                                      
    WGS_Epistasis.extend(expand('../02.mapping/mosdepth_coverage/{sample}.mosdepth.global.dist.txt',
                                           sample=samples.keys())) 
    WGS_Epistasis.extend(expand('../02.mapping/mosdepth_coverage/{sample}.mosdepth.summary.txt',
                                           sample=samples.keys())) 
    WGS_Epistasis.extend(expand('../02.mapping/qualimap_report/{sample}/qualimapReport.html',
                                           sample=samples.keys())) 
    WGS_Epistasis.extend(expand('../02.mapping/qualimap_report/{sample}/genome_results.txt',
                                           sample=samples.keys())) 
    WGS_Epistasis.extend(expand('../02.mapping/samtools_flagstat/{sample}_dup_bam_flagstat.tsv',
                                           sample=samples.keys())) 
    WGS_Epistasis.extend(expand('../02.mapping/samtools_stats/{sample}_dup_bam_stats.tsv',
                                           sample=samples.keys())) 
    # Call Variant
    WGS_Epistasis.extend(expand('../03.call_variant/{sample}.vcf.gz',
                                           sample=samples.keys()))
    if config['print_target']:
        logger.info(WGS_Epistasis)
    return  WGS_Epistasis
    
def judge_bwa_index(config:dict = None) -> bool:
    """
    判断是否需要重新构建bwa索引
    """
    bwa_index = config['bwa_mem2']['index']
    bwa_index_files = [bwa_index + suffix for suffix in ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac', '.alt']]
    
    return not all(os.path.exists(f) for f in bwa_index_files)

# --------------------- #
