use anyhow::{Context, Result};
use clap::Parser;
use regex::Regex;
// --- 【修改 1/4】: 引入 serde 用于 JSON 序列化 ---
use serde::Serialize;
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use walkdir::WalkDir;

/// 用于解析文件名的结构体，包含样本名和R1/R2信息
#[derive(Debug)]
struct SampleFileInfo {
    sample_name: String,
    read_pair: String,
    original_path: PathBuf,
}

// --- 【修改 2/4】: 新增用于生成 JSON 报告的结构体 ---
#[derive(Serialize, Debug)]
struct RenamingReportEntry {
    sample_name: String,
    // 使用相对于输出目录的路径
    new_r1_path_relative: String,
    new_r2_path_relative: String,
    // 使用绝对路径以方便溯源
    original_r1_path_absolute: String,
    original_r2_path_absolute: String,
    md5_r1: Option<String>,
    md5_r2: Option<String>,
}


/// CLI 参数定义
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(name = "seq_preprocessor")]
#[command(about = "自动整理不同来源的测序数据，统一命名和目录结构。")]
struct Cli {
    /// 原始数据所在的根目录路径
    #[arg(short, long)]
    input: PathBuf,

    /// 整理后数据的输出目录路径
    #[arg(short, long)]
    output: PathBuf,

    /// 指定在每个样本文件夹中生成的 MD5 文件的名称
    #[arg(long, default_value = "md5.txt")]
    md5_name: String,

    /// 在输出目录顶层生成一个包含所有文件信息的总 MD5 文件。
    /// 例如: --summary-md5 checksums.txt
    #[arg(long)]
    summary_md5: Option<PathBuf>,

    /// 一个开关，用于禁止在每个样本子目录中创建独立的 MD5 文件。
    #[arg(long)]
    no_per_sample_md5: bool,
    
    // --- 【修改 3/4】: 添加新的命令行参数 ---
    /// 生成一个 JSON 格式的重命名报告文件。
    #[arg(long)]
    json_report: Option<PathBuf>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if !cli.input.exists() {
        anyhow::bail!("输入路径不存在: {}", cli.input.display());
    }

    // 准备输出目录
    fs::create_dir_all(&cli.output).context(format!("无法创建输出目录: {}", cli.output.display()))?;

    println!("- 开始扫描输入目录: {}", cli.input.display());

    // 1. 收集所有 FASTQ 和 MD5 文件信息
    let mut fastq_files: Vec<SampleFileInfo> = Vec::new();
    let mut md5_files: Vec<PathBuf> = Vec::new();
    let mut unmatched_fastq_files: Vec<PathBuf> = Vec::new();

    // --- 【核心修改】: 使用两个独立的、更简单的正则表达式，提高鲁棒性 ---
    let re_illumina = Regex::new(r"^(.*?)_S\d+_L\d+_([Rr][12])_\d+\.f(ast)?q\.gz$").unwrap();
    let re_generic = Regex::new(r"^(.*)[\._]([Rr][12]|[12])\.f(ast)?q\.gz$").unwrap();


    for entry in WalkDir::new(&cli.input).into_iter().filter_map(|e| e.ok()) {
        let path = entry.path();
        if path.is_file() {
            let file_name = match path.file_name().and_then(|s| s.to_str()) {
                Some(name) => name,
                None => continue,
            };

            let mut matched = false;
            // --- 【核心修改】: 使用 if / else if 结构，按顺序尝试匹配 ---
            if let Some(caps) = re_illumina.captures(file_name) {
                if let (Some(sample), Some(pair_cap)) = (caps.get(1), caps.get(2)) {
                    let read_pair = if pair_cap.as_str().to_lowercase() == "r1" { "R1" } else { "R2" }.to_string();
                    fastq_files.push(SampleFileInfo {
                        sample_name: sample.as_str().to_string(),
                        read_pair,
                        original_path: path.to_path_buf(),
                    });
                    matched = true;
                }
            } else if let Some(caps) = re_generic.captures(file_name) {
                if let (Some(sample), Some(pair_cap)) = (caps.get(1), caps.get(2)) {
                    let read_pair = match pair_cap.as_str().to_lowercase().as_str() {
                        "r1" | "1" => "R1".to_string(),
                        "r2" | "2" => "R2".to_string(),
                        _ => unreachable!(),
                    };
                    fastq_files.push(SampleFileInfo {
                        sample_name: sample.as_str().to_string(),
                        read_pair,
                        original_path: path.to_path_buf(),
                    });
                    matched = true;
                }
            }
            
            if !matched {
                 if file_name.ends_with(".fq.gz") || file_name.ends_with(".fastq.gz") {
                    unmatched_fastq_files.push(path.to_path_buf());
                } else if file_name.to_lowercase().contains("md5") && file_name.ends_with(".txt") {
                    md5_files.push(path.to_path_buf());
                }
            }
        }
    }

    if !unmatched_fastq_files.is_empty() {
        let mut error_message =
            "错误：发现以下 FASTQ 文件命名不符合预期的两种模式，请检查或修正文件名:\n".to_string();
        for path in unmatched_fastq_files {
            error_message.push_str(&format!("  - {}\n", path.display()));
        }
        error_message.push_str("\n预期的模式为: \n  1. <样本名>_<Read号>.fq.gz (例如: B11_1.fq.gz)\n  2. <样本名>.<Read号>.fq.gz (例如: CK_L_10_R1.R1.fq.gz)\n  3. Illumina 格式 (例如: Sample_S1_L001_R1_001.fastq.gz)\n");
        anyhow::bail!(error_message);
    }

    println!("- 扫描完成: 找到 {} 个符合规范的 FASTQ 文件, {} 个 MD5 文件。", fastq_files.len(), md5_files.len());
    
    let mut checksum_map: HashMap<String, String> = HashMap::new();
    println!("- 开始解析 MD5 文件...");
    for md5_file_path in &md5_files {
        let file = fs::File::open(md5_file_path)?;
        let reader = BufReader::new(file);
        for line in reader.lines().filter_map(|l| l.ok()) {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 2 { continue; }
            let checksum = parts[0].to_string();
            let original_filename = PathBuf::from(parts[1])
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string();
            checksum_map.insert(original_filename, checksum);
        }
    }
    println!("- MD5 解析完成，共加载 {} 条记录。", checksum_map.len());


    let mut samples: HashMap<String, (Option<PathBuf>, Option<PathBuf>)> = HashMap::new();
    for file_info in fastq_files {
        let entry = samples.entry(file_info.sample_name.clone()).or_default();
        if file_info.read_pair == "R1" {
            entry.0 = Some(file_info.original_path);
        } else if file_info.read_pair == "R2" {
            entry.1 = Some(file_info.original_path);
        }
    }

    println!("- 已将文件整理为 {} 个独立样本。", samples.len());
    
    let mut summary_md5_lines: Vec<String> = Vec::new();
    // --- 【修改 4/4】: 初始化一个 Vec 来收集 JSON 报告的条目 ---
    let mut json_report_entries: Vec<RenamingReportEntry> = Vec::new();

    for (sample_name, (r1_opt, r2_opt)) in samples {
        println!("  - 正在处理样本: {}", sample_name);

        let (original_r1, original_r2) = match (r1_opt, r2_opt) {
            (Some(r1), Some(r2)) => (r1, r2),
            _ => {
                eprintln!("  ! 警告: 样本 {} 缺少 R1 或 R2 文件，已跳过。", sample_name);
                continue;
            }
        };

        let sample_output_dir = cli.output.join(&sample_name);
        fs::create_dir_all(&sample_output_dir)
            .context(format!("无法为样本 {} 创建目录", sample_name))?;

        let new_r1_name = format!("{}_R1.fq.gz", sample_name);
        let new_r2_name = format!("{}_R2.fq.gz", sample_name);
        let new_r1_path = sample_output_dir.join(&new_r1_name);
        let new_r2_path = sample_output_dir.join(&new_r2_name);

        #[cfg(unix)]
        {
            if new_r1_path.exists() { fs::remove_file(&new_r1_path)?; }
            if new_r2_path.exists() { fs::remove_file(&new_r2_path)?; }
            std::os::unix::fs::symlink(&original_r1, &new_r1_path)
                .context(format!("无法创建软链接: {}", new_r1_path.display()))?;
            std::os::unix::fs::symlink(&original_r2, &new_r2_path)
                .context(format!("无法创建软链接: {}", new_r2_path.display()))?;
        }
        #[cfg(not(unix))]
        {
             fs::copy(&original_r1, &new_r1_path)
                .context(format!("无法复制文件: {}", new_r1_path.display()))?;
             fs::copy(&original_r2, &new_r2_path)
                .context(format!("无法复制文件: {}", new_r2_path.display()))?;
        }

        let original_r1_filename = original_r1.file_name().unwrap().to_str().unwrap();
        let original_r2_filename = original_r2.file_name().unwrap().to_str().unwrap();
        
        let checksum_r1 = checksum_map.get(original_r1_filename).cloned();
        let checksum_r2 = checksum_map.get(original_r2_filename).cloned();

        if !cli.no_per_sample_md5 {
            let mut per_sample_md5_content = String::new();
            if let Some(c) = &checksum_r1 { per_sample_md5_content.push_str(&format!("{}  {}\n", c, new_r1_name)); }
            if let Some(c) = &checksum_r2 { per_sample_md5_content.push_str(&format!("{}  {}\n", c, new_r2_name)); }

            if !per_sample_md5_content.is_empty() {
                let per_sample_md5_path = sample_output_dir.join(&cli.md5_name);
                fs::write(&per_sample_md5_path, per_sample_md5_content)
                    .context(format!("无法写入样本 MD5 文件: {}", per_sample_md5_path.display()))?;
                println!("    - 已生成样本 MD5 文件: {}", cli.md5_name);
            }
        }
        
        let relative_r1_path = PathBuf::from(&sample_name).join(&new_r1_name);
        let relative_r2_path = PathBuf::from(&sample_name).join(&new_r2_name);
        if let Some(c) = &checksum_r1 { summary_md5_lines.push(format!("{}  {}", c, relative_r1_path.display())); }
        if let Some(c) = &checksum_r2 { summary_md5_lines.push(format!("{}  {}", c, relative_r2_path.display())); }
        
        // --- 填充 JSON 报告条目 ---
        if cli.json_report.is_some() {
            json_report_entries.push(RenamingReportEntry {
                sample_name: sample_name.clone(),
                new_r1_path_relative: relative_r1_path.to_string_lossy().to_string(),
                new_r2_path_relative: relative_r2_path.to_string_lossy().to_string(),
                original_r1_path_absolute: fs::canonicalize(&original_r1)?.to_string_lossy().to_string(),
                original_r2_path_absolute: fs::canonicalize(&original_r2)?.to_string_lossy().to_string(),
                md5_r1: checksum_r1,
                md5_r2: checksum_r2,
            });
        }
    }
    
    if let Some(summary_md5_filename) = &cli.summary_md5 {
        let summary_path = cli.output.join(summary_md5_filename);
        if !summary_md5_lines.is_empty() {
            summary_md5_lines.sort();
            let final_content = summary_md5_lines.join("\n") + "\n";
            fs::write(&summary_path, final_content)
                .context(format!("无法写入总纲 MD5 文件: {}", summary_path.display()))?;
            println!("\n- 成功生成总纲 MD5 文件于: {}", summary_path.display());
        } else {
            println!("\n- 未找到任何 MD5 信息，因此未生成总纲 MD5 文件。");
        }
    }
    
    // --- 在所有样本处理完后，写入 JSON 报告文件 ---
    if let Some(report_path) = &cli.json_report {
        if !json_report_entries.is_empty() {
            println!("\n- 正在生成 JSON 报告文件...");
            // 使用 serde_json 将数据结构转换为格式化的 JSON 字符串
            let report_json = serde_json::to_string_pretty(&json_report_entries)?;
            fs::write(report_path, report_json)
                .context(format!("无法写入 JSON 报告文件: {}", report_path.display()))?;
            println!("- 成功生成 JSON 报告文件于: {}", report_path.display());
        } else {
             println!("\n- 未找到任何样本，因此未生成 JSON 报告文件。");
        }
    }

    println!("\n- 所有任务完成！标准化的数据已存放于: {}", cli.output.display());
    println!("- 提示: 脚本默认使用软链接（symlink）来指向原始文件，这不会消耗额外磁盘空间。");

    Ok(())
}
