***Author  : JZHANG***  
***Date    : 2025-10-17***  
***Version : 1.0v***
# Epistasis Analysis
## Introduction


## Workflow
### Quality control
`md5_checker_rs`
为了加快分析速度在`scripts/md5_checker_rs`文件夹中使用`rust`语言编写了下机数据`MD5`质控脚本。`target`文件夹下的可直接运行的二进制文件为`md5_checker_rs`。是在`Ubuntu 20.04.4 LTS x86_64`系统下编译，如果需要在不同的平台下使用，请自行编译,编译命令:`cargo build --release`。
```bash
./md5_checker_rs --help
使用多个线程自动查找和验证指定目录中的 MD5 校验和文件。

Usage: md5_checker_rs [OPTIONS] <DIRECTORIES>...

Arguments:
  <DIRECTORIES>...  一个或多个包含数据文件和 MD5 校验和文件的目录路径。

Options:
  -f, --filename <FILENAME>  MD5 校验和文件的名称。 [default: MD5.txt]
  -t, --threads <THREADS>    用于并发验证的线程数 (0 表示使用 Rayon 的默认值，通常是 CPU 核心数)。 [default: 0]
  -o, --output <OUTPUT>      生成 TSV 格式报告文件的路径。
  -h, --help                 Print help
  -V, --version              Print version

示例: md5_checker_rs /path/to/data1 /path/to/data2 -f checksums.txt -t 16 -o report.tsv
```
### Run snakemake pipeline
```bash
# Creating conda environments via mamba
snakemake --use-conda --conda-create-envs-only --conda-frontend mamba
# Clean install conda/mamba envs
snakemake --use-conda --conda-cleanup-envs
# Runing the snakemake pipeline via mamba
snakemake --cores=50 -p --conda-frontend mamba --use-conda
```