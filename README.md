***Author  : JZHANG***  
***Date    : 2025-11-1***  
***Version : 2.0v***
# Epistasis Analysis
## Introduction
`Epistasis Analysis`项目是基于WGS+RNA数据对群体中存在的上位性进行研究的代码仓库
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
## Run snakemake pipeline

### Install `snakemake=9.9.0`
```bash
# fix conda channel_priority -> flexible
conda config --set channel_priority flexible
# Create snakemake conda environment and install snakemak
conda create -n snakemake
# activate snakemake
conda activate snakemake
# install snakemake=9.9.0
conda install bioconda::snakemake=9.9.0
```
### Deployment Snakemake analysis environment via mamba & run pipeline
```bash
# Creating conda environments via mamba
snakemake --use-conda --conda-create-envs-only --conda-frontend mamba
# Clean install conda/mamba envs
snakemake --use-conda --conda-cleanup-envs
# Runing the snakemake pipeline via mamba
snakemake --cores=50 -p --conda-frontend mamba --use-conda
```
```bash
# 测试 & 开发模式
snakemake --dry-run
```

### 运行 Nextflow 流程 (Run Nextflow pipeline)
本流程同样支持 Nextflow (DSL2) 进行部署和运行。
#### 1.安装 Nextflow
首先，你需要安装 Nextflow。根据 nextflow.config 中的 manifest 定义，推荐使用 v23.04.0 或更高版本。
```bash
# 推荐使用 SDKMAN 安装：
curl -s "https://get.sdkman.io" | bash
# (打开新终端后)
sdk install nextflow
# 或者使用 cURL 直接下载：
curl -s https://get.nextflow.io | bash
# (将 nextflow 可执行文件移动到 $PATH 路径下)
```
#### 2. 理解 profiles 配置文件
本流程使用 profiles 来管理执行器（在哪里运行）和软件环境（如何运行）。你可以在 conf/profiles.config 中查看和修改它们的定义。
主要 profile 分类:
  执行器 (Executor):
  - local: 在你的本地机器上执行流程（默认并发数较高，请见配置）。
  - slurm: 将任务提交到 Slurm 集群调度系统。
  软件管理器 (Software):
  - conda: 使用 Conda / Micromamba 来自动创建和管理 envs/ 目录下的环境。
  - docker: 使用 Docker 容器来运行（需要 Docker 服务运行中）。
  - singularity: slurm profile 默认启用了 Singularity (Apptainer) 容器支持。
  你必须组合使用这些 profile 来指定运行方式。

#### 3.运行流程
基础运行 (Conda)
这是你提供的基础运行命令，使用 conda profile，并连接到 Nextflow Tower 进行监控：
```bash
# 注意: 这个命令依赖于 nextflow.config 中设置的默认执行器
./Epistasis.nf -resume -with-tower -name "epistasis_$(date +%Y%m%d_%H%M%S)" -profile conda
```

推荐的运行方式 (组合 Profile)
为了更明确地控制执行器和环境，推荐组合使用 profile，并使用 --参数 来指定输入/输出：


```bash
# 本地运行 (使用 Conda):
./Epistasis.nf -resume -profile local,conda --outdir ./results --pop "EUR_cohort"
# 集群运行 (使用 Singularity - Slurm profile 内置):
./Epistasis.nf -resume -profile slurm --outdir /path/to/results --pop "ASN_cohort"
# 本地运行 (使用 Docker):
./Epistasis.nf -resume -profile local,docker --outdir ./results --pop "AFR_cohort"
```
#### 4.主要参数 (Main Parameters)
你可以在命令行中使用 --参数名 <值> 的方式来指定 Epistasis.nf 脚本中的 params。
  --input: (必需) 输入样本表，例如 samplesheet.csv。
  --outdir: (可选) 所有输出结果的保存目录。 [默认: ./results]
  --pop: (可选) 群体的名称，用于命名输出文件或目录。 [默认: default_population]
## Reference