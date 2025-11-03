# Epistasis Analysis Pipeline

本项目是一个为上位性（Epistasis）分析设计的高性能全基因组测序（WGS）数据处理流程。它内置了从原始测序数据（FASTQ）到经过质控、过滤和注释的变异位点（VCF）的完整分析步骤。

为了提供最大的灵活性和可复现性，该流程同时使用了 **Nextflow** 和 **Snakemake** 两种主流生物信息工作流语言进行实现。

此外，项目还包含两个使用 Rust 开发的命令行工具 `seq_preprocessor` 和 `json_md5_verifier`，用于在分析开始前实现高效、安全的数据整理和校验。

## 核心特性

- **双引擎驱动**：提供 Nextflow (`Epistasis.nf`) 和 Snakemake (`snakefile`) 两套完全独立的流程实现，用户可根据自己的熟悉程度和计算环境任选其一。
- **标准化分析步骤**：流程集成了当前主流的生物信息学分析工具，包括 `Fastp`, `BWA-MEM2`, `Samtools`, `BCFtools`, `SnpEff` 等。
- **高效并行计算**：利用 `bcftools` 按染色体拆分任务，并行进行变异检测，显著加快分析速度。
- **环境隔离与可复现**：通过 Conda (`Snakemake`) 和 Micromamba (`Nextflow`) 自动管理软件环境，确保流程在不同系统中的可复现性。
- **智能数据预处理**: 内置 Rust 工具，用于自动化整理和校验原始数据。该步骤在 Snakemake 中全自动执行，在 Nextflow 中可手动调用。

## 目录结构

```
.
├── config/               # Nextflow 和 Snakemake 的详细配置文件
├── data/                 # 存放通用数据，如 Adapter 序列
├── envs/                 # Snakemake 和 Nextflow 共用的 Conda 环境定义 (YAML)
├── modules/              # Nextflow 流程的模块脚本 (.nf)
├── rules/                # Snakemake 流程的规则脚本 (.smk)
├── scripts/              # 存放定制开发的辅助脚本，包括 Rust 工具源码
├── docs/                 # 存放多语言文档
│   └── README_CN.md      # 中文文档
├── Epistasis.nf          # Nextflow 主流程文件
├── nextflow.config       # Nextflow 主配置文件
├── snakefile             # Snakemake 主流程文件
├── config.yaml           # Snakemake 主配置文件
├── sample.csv            # Snakemake 使用的样本表示例文件
├── @sample_nextflow.csv  # Nextflow 使用的样本表示例文件
└── README.md             # 英文文档 (English documentation)
```

---

## 工作流程

根据您选择的流程引擎（Nextflow 或 Snakemake），数据准备的步骤有所不同。

- **对于 Snakemake**：数据整理和校验步骤已**完全内嵌**到工作流中，实现了从原始数据到最终结果的“一键运行”。
- **对于 Nextflow**：需要**手动执行**数据预处理步骤，以生成符合流程要求的标准化输入文件。

---

### **自定义数据处理工具**

本项目提供了两个使用 Rust 开发的高性能命令行工具，用于数据预处理。对于 Snakemake 流程，这些工具会被自动调用；对于 Nextflow 流程，您需要手动运行它们。

#### 1. 编译 Rust 工具 (首次使用)

```bash
# 编译 seq_preprocessor
cd scripts/seq_preprocessor && cargo build --release

# 编译 json_md5_verifier
cd ../json_md5_verifier && cargo build --release
```
编译后的可执行文件位于各自的 `target/release/` 目录下。建议将它们复制到您的 `PATH` 路径下以便于调用。

#### 2. 工具介绍

- **`seq_preprocessor`**: 此工具能扫描原始数据目录，自动识别并整理 FASTQ 文件，最终通过软链接（symlink）的方式创建一个标准化的样本目录，同时生成一份包含文件路径和 MD5 信息的 JSON 报告。

- **`json_md5_verifier`**: 此工具读取上述 JSON 报告，并多线程并发校验文件的 MD5 值，以确保数据在整理过程中的完整性。

---

### **运行分析流程**

#### 1. 依赖安装

- **核心流程引擎**:
    - **Nextflow**: [安装 Java 8 或更高版本](https://www.nextflow.io/docs/latest/getstarted.html#requirements)，然后运行 `curl -s https://get.nextflow.io | bash`。
    - **Snakemake**: `conda install -c conda-forge -c bioconda snakemake`
- **环境管理工具**:
    - **Snakemake 用户**: 需要安装 `conda` 或 `mamba`。
    - **Nextflow 用户**: 需要安装 `micromamba`。

#### 2. 参数配置与运行

##### **使用 Snakemake (自动化数据准备)**

1.  **配置**: 编辑 `config.yaml` 文件。您**只需**在 `raw_data_path` 下指定一个或多个存放**原始 FASTQ 文件**的目录路径。同时，配置好 `genome` (参考基因组)等其他关键路径。

2.  **运行**: Snakemake 将首先自动调用 `seq_preprocessor` 和 `json_md5_verifier` 处理您的原始数据，然后继续执行后续的 WGS 分析。

    ```bash
    # 运行流程 (确保已安装 conda/mamba)
    snakemake --snakefile snakefile \
              --cores 32 \
              --use-conda
    ```

##### **使用 Nextflow (手动数据准备)**

1.  **手动准备数据**: 
    首先，使用 `seq_preprocessor` 和 `json_md5_verifier` 工具处理您的原始数据。假设整理后的数据存放在 `./wgs_input` 目录。

    ```bash
    # 步骤1: 整理数据
    seq_preprocessor --input /path/to/raw_data --output ./wgs_input --json-report ./wgs_input/report.json

    # 步骤2: 校验数据
    json_md5_verifier --input ./wgs_input/report.json --base-dir ./wgs_input
    ```

2.  **创建样本文件**: 参照 `@sample_nextflow.csv` 的格式，创建一个样本文件，路径指向**整理后**的文件。

3.  **配置**: 修改 `nextflow.config` 中的 `params` 部分，特别是 `bcftools_reference` 等路径。

4.  **运行**:
    ```bash
    # 运行流程 (确保已安装 micromamba)
    nextflow run Epistasis.nf \
        -profile conda \
        --input your_sample_sheet.csv \
        --outdir ./results_nextflow
    ```

- `--cores` / `-profile conda`: 指定核心数或使用 conda/micromamba 环境。
- `--input` / `--outdir`: 指定输入样本表和输出目录。

##### **使用 Nextflow Tower (可选)**

本流程支持使用 [Nextflow Tower](https://tower.nf/) 进行监控和管理。

1.  **获取 Token**: 访问 [tower.nf](https://tower.nf/) 并登录，在 "Your Tokens" 部分创建一个新的 Token。
2.  **配置环境变量**: 在您的终端环境中，导出 Tower Token 和 Endpoint：
    ```bash
    export TOWER_ACCESS_TOKEN='YOUR_TOKEN_HERE'
    export TOWER_ENDPOINT='https://api.tower.nf'
    ```
3.  **运行**: 在您的 `nextflow run` 命令中加入 `-with-tower` 参数。
    ```bash
    nextflow run Epistasis.nf \
        -with-tower \
        -profile conda \
        ...
    ```
    现在，您可以在 Tower 网站上实时监控流程的执行情况。

## 流程产出

分析完成后，您将在指定的输出目录（例如 `results_nextflow`）中找到一系列结构化的子目录，包含各个分析步骤的结果：

- `01.fastp_clean/`: 质控后的 FASTQ 文件和报告。
- `02.mapping/`: BAM 比对文件、索引及各种 mapping 质控报告。
- `03.variant_calling/`: VCF 文件（包括原始、合并、过滤和注释后的版本）及相关统计。
- `multiqc/`: 聚合所有分析步骤的 MultiQC 综合报告。

您可以从 `multiqc/` 中的 HTML 报告开始，全面评估数据质量。最终，`03.variant_calling/` 中经过滤和注释的 VCF 文件可用于下游的上位性分析。
