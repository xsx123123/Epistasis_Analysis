#!/usr/bin/env rust-script
// ... (use 语句保持不变) ...
use std::fs::{self, File};
use std::io::{self, BufReader, Read};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::env;

use anyhow::Result;
use clap::Parser;
use log::{error, info, warn, LevelFilter};
use rayon::prelude::*;

use indicatif::{ProgressBar, ProgressStyle, ParallelProgressIterator};

/// 使用多个线程自动查找和验证指定目录中的 MD5 校验和文件。
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None, after_help = "示例: md5_checker_rs /path/to/data1 /path/to/data2 -f checksums.txt -t 16 -o report.tsv --log-file logs/checker.log")]
struct Cli {
    /// 一个或多个包含数据文件和 MD5 校验和文件的目录路径。
    #[arg(required = true)]
    directories: Vec<PathBuf>,

    /// MD5 校验和文件的名称。
    #[arg(short, long, default_value = "MD5.txt")]
    filename: String,

    /// 用于并发验证的线程数 (0 表示使用 Rayon 的默认值，通常是 CPU 核心数)。
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// 生成 TSV 格式报告文件的路径。
    #[arg(short, long)]
    output: Option<PathBuf>,

    // ----- 【修改 1/3】: 在这里添加新的命令行参数 -----
    /// 指定日志文件的完整路径。
    #[arg(long, default_value = "md5_checker.log")]
    log_file: PathBuf,
}

// ... (VerificationTask 和 VerificationResult 结构体不变) ...
#[derive(Debug, Clone)]
struct VerificationTask {
    expected_md5: String,
    file_to_check: PathBuf,
    relative_filename: String,
}

#[derive(Debug)]
struct VerificationResult {
    timestamp: String,
    full_path: String,
    expected_md5: String,
    actual_md5: String,
    status: &'static str,
    message: String,
}


// ----- 【修改 2/3】: 修改函数签名，使其接收一个路径参数 -----
/// 配置日志记录器，同时输出到控制台和文件。
fn setup_logger(log_path: &Path) -> Result<()> {
    // 确保日志文件所在的目录存在
    if let Some(parent_dir) = log_path.parent() {
        fs::create_dir_all(parent_dir)?;
    }

    fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "{} - {} - {}",
                chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
                record.level(),
                message
            ))
        })
        .level(LevelFilter::Info)
        .chain(io::stdout())
        // ----- 使用传入的参数作为日志文件名 -----
        .chain(fern::log_file(log_path)?)
        .apply()?;
    Ok(())
}

// ... (calculate_md5, verify_file_task, generate_report 函数不变) ...
fn calculate_md5(filepath: &Path) -> io::Result<String> {
    let file = File::open(filepath)?;
    let mut reader = BufReader::new(file);
    let mut context = md5::Context::new();
    let mut buffer = [0; 65536]; 

    loop {
        let bytes_read = reader.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        context.consume(&buffer[..bytes_read]);
    }

    Ok(format!("{:x}", context.compute()))
}

fn verify_file_task(task: &VerificationTask) -> VerificationResult {
    let timestamp = chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
    let full_path_str = task.file_to_check.to_string_lossy().to_string();

    if !task.file_to_check.exists() {
        return VerificationResult {
            timestamp,
            full_path: full_path_str,
            expected_md5: task.expected_md5.clone(),
            actual_md5: "N/A".to_string(),
            status: "FAIL",
            message: "File not found".to_string(),
        };
    }

    match calculate_md5(&task.file_to_check) {
        Ok(actual_md5) => {
            if actual_md5 == task.expected_md5 {
                VerificationResult { timestamp, full_path: full_path_str, expected_md5: task.expected_md5.clone(), actual_md5, status: "PASS", message: "MD5 match".to_string() }
            } else {
                VerificationResult { timestamp, full_path: full_path_str, expected_md5: task.expected_md5.clone(), actual_md5, status: "FAIL", message: "MD5 mismatch".to_string() }
            }
        }
        Err(e) => VerificationResult { timestamp, full_path: full_path_str, expected_md5: task.expected_md5.clone(), actual_md5: "N/A".to_string(), status: "FAIL", message: format!("Could not calculate MD5 (read error: {})", e) },
    }
}

fn generate_report(results: &[VerificationResult], output_file: &Path) -> Result<()> {
    info!("--- 正在生成验证报告至: {} ---", output_file.display());
    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(output_file)?;
    writer.write_record(&["CheckTime", "FilePath", "ExpectedMD5", "ActualMD5", "Status", "Message"])?;
    for res in results {
        writer.write_record(&[&res.timestamp, &res.full_path, &res.expected_md5, &res.actual_md5, res.status, &res.message])?;
    }
    writer.flush()?;
    info!("报告生成完成。");
    Ok(())
}


fn main() -> Result<()> {
    let cli = Cli::parse();
    
    // ----- 【修改 3/3】: 将解析出的路径传递给 setup_logger -----
    setup_logger(&cli.log_file)?;

    if cli.threads > 0 {
        unsafe {
        env::set_var("RAYON_NUM_THREADS", cli.threads.to_string());
        }
    }

    // ... (main 函数的其余部分保持不变) ...
    let mut tasks = Vec::new();
    info!("--- 开始搜索 {} 文件并构建任务列表 ---", cli.filename);
    for directory in &cli.directories {
        if !directory.is_dir() {
            warn!("提供的路径不是一个目录，已跳过: {}", directory.display());
            continue;
        }
        let md5_file_path = directory.join(&cli.filename);
        if !md5_file_path.exists() {
            warn!("在目录 {} 中未找到 {}，已跳过。", directory.display(), cli.filename);
            continue;
        }
        info!("找到并正在解析: {}", md5_file_path.display());
        let content = fs::read_to_string(&md5_file_path)?;
        for (i, line) in content.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() { continue; }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 2 {
                warn!("警告: 在 {} 中发现格式错误的行 #{}, 已跳过: '{}'", md5_file_path.display(), i + 1, line);
                continue;
            }
            let expected_md5 = parts[0].to_lowercase();
            let relative_filename = parts[1].to_string();
            let file_to_check = directory.join(&relative_filename);
            tasks.push(VerificationTask { expected_md5, file_to_check, relative_filename });
        }
    }

    let num_tasks = tasks.len();
    if num_tasks == 0 {
        warn!("未找到需要验证的文件。脚本执行完毕。");
        return Ok(());
    }

    info!(
        "--- 找到 {} 个文件待验证，开始使用 {} 个线程进行处理 ---",
        num_tasks,
        rayon::current_num_threads()
    );

    let pb = ProgressBar::new(num_tasks as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) - {per_sec} - ETA: {eta}")?
        .progress_chars("#>-"));

    let has_failures = AtomicBool::new(false);
    
    let results: Vec<VerificationResult> = tasks
        .par_iter()
        .progress_with(pb)
        .map(|task| {
            let result = verify_file_task(task);
            if result.status == "FAIL" {
                has_failures.store(true, Ordering::Relaxed);
                error!("  [FAIL] {}: {}", result.message, task.relative_filename);
                if result.message == "MD5 mismatch" {
                    error!("      预期: {}", result.expected_md5);
                    error!("      实际:   {}", result.actual_md5);
                }
            } 
            result
        })
        .collect();

    if let Some(output_path) = cli.output {
        if let Err(e) = generate_report(&results, &output_path) {
            error!("无法写入报告文件 {}: {}", output_path.display(), e);
        }
    }

    info!("======================================================");
    if has_failures.load(Ordering::Relaxed) {
        error!("验证过程中发现错误。请检查日志和报告文件以获取详细信息。");
        info!("======================================================");
        std::process::exit(1);
    } else {
        info!("所有 {} 个文件均成功通过验证！", results.len());
        info!("======================================================");
    }

    Ok(())
}