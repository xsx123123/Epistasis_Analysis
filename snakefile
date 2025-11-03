#!/usr/bin/env python3
# *---utf-8---*
# Version: 0.1.1.a
# Author : JZHANG
# ------- snakemake version check ------- #
from snakemake.utils import min_version
min_version("9.9.0")
# --------- main snakefile --------- #
configfile: "config/config.yaml"
configfile: "config.yaml"
# --------- snakemake rule --------- #
# include all rules from the rules directory
include: 'rules/wgs/00.log.smk'
include: 'rules/wgs/00.common.smk'
include: 'rules/wgs/00.id_convert.smk'
include: 'rules/wgs/01.file_convert_md5.smk'
include: 'rules/wgs/01.short_read_qc.smk'
include: 'rules/wgs/01.Contamination_check.smk'
include: 'rules/wgs/01.short_read_clean.smk'
include: 'rules/wgs/02.mapping.smk'
include: 'rules/wgs/03.callvraiant_help.smk'
include: 'rules/wgs/03.callvraiant.smk'
include: 'rules/wgs/F.multi_merge_qc.smk'
# --------- target rule --------- #
rule all:
    input:
        Epistasis(config = config)
# --------- target rule --------- #