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
    short-read & long-read RNA-seq assembly & annotation pipeline target
    """
    # short-read raw-data qc result
    WGS_Epistasis = [
            "../01.qc/md5_check.tsv",
            "../01.qc/short_read_r1_multiqc/",
            "../01.qc/short_read_r2_multiqc/"    
                    
        ]
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

    WGS_Epistasis.extend(expand("../02.mapping/bwa_mem2/{sample}.sort.bam",
                                           sample=samples.keys()))
    WGS_Epistasis.extend(expand("../02.mapping/bwa_mem2/{sample}.sort.bam.bai",
                                           sample=samples.keys()))
    WGS_Epistasis.extend(expand('../02.mapping/bwa_mem2/{sample}.dup.sort.bam',
                                           sample=samples.keys()))
    WGS_Epistasis.extend(expand('../02.mapping/bwa_mem2/{sample}_marked_dup_metrics.txt',
                                           sample=samples.keys()))
    WGS_Epistasis.extend(expand('../02.mapping/bwa_mem2/{sample}.dup.sort.bam.bai',
                                           sample=samples.keys()))                                      
    WGS_Epistasis.extend(expand('../02.mapping/mosdepth_coverage/{sample}.mosdepth.global.dist.tx',
                                           sample=samples.keys())) 
    WGS_Epistasis.extend(expand('../02.mapping/mosdepth_coverage/{sample}.mosdepth.summary.txt',
                                           sample=samples.keys())) 
    WGS_Epistasis.extend(expand('../02.mapping/qualimap_report/{sample}_qualimap_report.html',
                                           sample=samples.keys())) 
    # Print Target rule           
    if config['print_target']:
        rich_print(WGS_Epistasis)
    return  WGS_Epistasis
    
# --------------------- #