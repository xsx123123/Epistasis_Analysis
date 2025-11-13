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
    Epistasis analysis function.
    This function performs epistasis
    analysis on the input configuration
    and returns a list of results.
    """
    # short-read raw-data qc result
    WGS_Epistasis = ["../01.qc/md5_check.tsv",
                    os.path.join('../00.raw_data',config['convert_md5']),
                    os.path.join('../00.raw_data',config['convert_md5'],"raw_data_md5.json")
                    ]
    # fastq-screen
    if config['fastq_screen']['run']:
        WGS_Epistasis.extend(expand("../01.qc/fastq_screen_r1/{sample}_R1_screen.txt",
                                          sample=samples.keys()))
        WGS_Epistasis.extend(expand("../01.qc/fastq_screen_r2/{sample}_R2_screen.txt",
                                          sample=samples.keys()))
        WGS_Epistasis.append("../01.qc/fastq_screen_multiqc_r1/multiqc_r1_fastq_screen_report.html")
        WGS_Epistasis.append("../01.qc/fastq_screen_multiqc_r2/multiqc_r2_fastq_screen_report.html")
    # fastqc & multiqc
    WGS_Epistasis.append("../01.qc/short_read_r1_multiqc/multiqc_r1_raw-data_report.html")
    WGS_Epistasis.append("../01.qc/short_read_r2_multiqc/multiqc_r2_raw-data_report.html")        
    # short-read trim & clean result
    WGS_Epistasis.append("../01.qc/multiqc_short_read_trim/multiqc_short_read_trim_report.html")
    # mapping result
    WGS_Epistasis.extend(expand("../02.mapping/bwa_mem2/temp/{sample}.raw.bam",
                                           sample=samples.keys()))
    WGS_Epistasis.extend(expand("../02.mapping/bwa_mem2/sort_index/{sample}.sort.bam",
                                           sample=samples.keys()))
    WGS_Epistasis.extend(expand("../02.mapping/bwa_mem2/MarkDup/{sample}.sort.Dup.bam",
                                           sample=samples.keys()))
    WGS_Epistasis.extend(expand("../02.mapping/bwa_mem2/MarkDup/{sample}.sort.Dup.bam.bai",
                                           sample=samples.keys()))
    WGS_Epistasis.extend(expand("../02.mapping/qualimap_report/{sample}/qualimapReport.html",
                                           sample=samples.keys())) 
    # WGS_Epistasis.extend(expand("../02.mapping/mosdepth_coverage/{sample}.sort.Dup.global.dist.txt",
    #                                        sample=samples.keys()))
    # WGS_Epistasis.extend(expand("../02.mapping/mosdepth_coverage/{sample}.sort.Dup.mosdepth.summary.txt",
    #                                        sample=samples.keys()))
    WGS_Epistasis.extend(expand("../02.mapping/samtools_stats/{sample}_dup_bam_stats.tsv",
                                           sample=samples.keys()))
    # Call Variant
    WGS_Epistasis.extend(expand('../03.call_variant/{sample}.raw.vcf.gz',
                                 sample=samples.keys()))
    WGS_Epistasis.append("../03.call_variant/variant_stats_multiqc/")
    WGS_Epistasis.append('../03.call_variant/merge.vcf.gz')
    WGS_Epistasis.append('../03.call_variant/merge_filter.sort.vcf.gz')
    # annotation vcf
    WGS_Epistasis.append('../03.call_variant/merge_filter.sort.annotation.csv')
    WGS_Epistasis.append('../03.call_variant/merge_filter.sort.annotation.html')
    WGS_Epistasis.append('../03.call_variant/merge_filter.sort.annotation.vcf') 
    # Report
    WGS_Epistasis.append('../F.report/Data_QC_report/Data_QC_report.html') 
    WGS_Epistasis.append('../F.report/Mapping_QC_report/Mapping_QC_report.html') 
    WGS_Epistasis.append("../F.report/Mapping_QC_report/Variant_Stats_report.html") 
    
    if config['print_target']:
       rich_print(WGS_Epistasis)
    return  WGS_Epistasis

def get_sample_data_dir(sample_id:str = None,
                        config:dict = None) -> str:
    """
    根据 *具体的* sample_id (e.g., "Sample_A"),
    查找 *包含* fastq 文件的 *目录*。
    
    注意：我修改了它，使其不再依赖 wildcards，
    而是直接接收 sample_id 字符串。
    """
    
    for base_dir in config["raw_data_path"]:
        sample_dir = os.path.join(base_dir, sample_id)
        if os.path.isdir(sample_dir):
            return sample_dir
                
    raise FileNotFoundError(f"无法在 {config['raw_data_path']} 中找到 {sample_id} 的数据目录")

def get_all_input_dirs(sample_keys:str = None,
                       config:dict = config) -> list:
    """
    遍历所有样本 ID，调用 get_sample_data_dir，
    返回一个包含所有数据目录的列表。
    """
    dir_list = []
    for sample_id in sample_keys:
        dir_list.append(get_sample_data_dir(sample_id,config = config))
    print(dir_list)
    return list(set(dir_list))

def judge_bwa_index(config:dict = None) -> bool:
    """
    判断是否需要重新构建bwa索引
    """
    bwa_index = config['bwa_mem2']['index']
    bwa_index_files = [bwa_index + suffix for suffix in ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac', '.alt']]
    
    return not all(os.path.exists(f) for f in bwa_index_files)

# --------------------- #
