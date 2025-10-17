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
    hybrid_rna_assembly = [
            "../01.qc/short_read_r1_multiqc/",
            "../01.qc/short_read_r2_multiqc/"            
        ]
    # short-read trim & clean result
    hybrid_rna_assembly.extend(expand("../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz",
                                          sample=load_samples.keys()))
    hybrid_rna_assembly.extend(expand("../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz",
                                          sample=load_samples.keys()))
    hybrid_rna_assembly.extend(expand("../01.qc/short_read_trim/{sample}.fastp.html",
                                          sample=load_samples.keys()))
    hybrid_rna_assembly.extend(expand("../01.qc/short_read_trim/{sample}.fastp.json",
                                          sample=load_samples.keys()))
    hybrid_rna_assembly.append("../01.qc/multiqc_short_read_trim/")
    # Print Target rule           
    if config['print_target']:
        rich_print(hybrid_rna_assembly)
    return  hybrid_rna_assembly


def _check_file_existence(path_str: str, name: str) -> None:
    """
    help function to check file existence
    """
    path = Path(path_str)
    if not path.exists():
        error_msg = f'Place check the {name} path: {path_str}'
        logger.error(error_msg)
        raise FileNotFoundError(error_msg)

def judge_file_optimized(config: Dict = None) -> None:
    """
    judge if the file is exist for config yaml file
    """
    if config is None:
        return

    software_paths = [
        (config["software"].get("tmhmm"), "tmhmm"),
        (config["software"].get("signalp"), "signalp"),
        (config["software"].get("interproscan"), "interproscan"),
    ]

    database_paths = [
        (config["busco"].get("lineage"), "BUSCO lineage database"),
        (config["mmseqs"].get("swissprot_database"), "swissprot database"),
        (config["annotate"].get("uniport_database"), "uniport database"),
        (config["annotate"].get("Pfam_A"), "Pfam A database"),
    ]

    all_paths = software_paths + database_paths

    for path_str, name in all_paths:
        if path_str:
            _check_file_existence(path_str, name)
