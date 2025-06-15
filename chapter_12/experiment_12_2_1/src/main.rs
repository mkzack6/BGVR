use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::collections::HashMap;

use anyhow::{Context, Result};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn, error};
use rayon::prelude::*;
use tokio::task;

/// Genomic Sparse Index Builder and Query Tool
#[derive(Parser, Debug)]
#[command(name = "rust_sparse_index")]
#[command(about = "A tool for building sparse genomic indices and querying allele frequencies")]
struct Args {
    /// Input genotype file path
    #[arg(short, long)]
    input: String,
    
    /// Output results file path
    #[arg(short, long, default_value = "index_results.txt")]
    output: String,
    
    /// Maximum genomic position index
    #[arg(short, long, default_value = "10000000")]
    max_index: usize,
    
    /// Number of parallel query tasks
    #[arg(short, long, default_value = "8")]
    num_tasks: usize,
    
    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,
}

/// Binary Indexed Tree (Fenwick Tree) for efficient range queries
#[derive(Debug, Clone)]
struct FenwickTree {
    data: Vec<u32>,
    size: usize,
}

impl FenwickTree {
    fn new(size: usize) -> Self {
        FenwickTree {
            data: vec![0; size + 1],
            size,
        }
    }

    fn update(&mut self, mut idx: usize, val: u32) -> Result<()> {
        if idx == 0 || idx > self.size {
            return Err(anyhow::anyhow!("Index {} out of bounds [1, {}]", idx, self.size));
        }
        
        while idx <= self.size {
            self.data[idx] = self.data[idx].saturating_add(val);
            idx += idx & (!idx + 1);
        }
        Ok(())
    }

    fn prefix_sum(&self, mut idx: usize) -> u32 {
        if idx > self.size {
            idx = self.size;
        }
        
        let mut result = 0u32;
        while idx > 0 {
            result = result.saturating_add(self.data[idx]);
            idx -= idx & (!idx + 1);
        }
        result
    }

    fn range_sum(&self, left: usize, right: usize) -> Result<u32> {
        if left > right || left == 0 || right > self.size {
            return Err(anyhow::anyhow!("Invalid range [{}, {}] for size {}", left, right, self.size));
        }
        
        let right_sum = self.prefix_sum(right);
        let left_sum = if left > 1 { self.prefix_sum(left - 1) } else { 0 };
        
        Ok(right_sum - left_sum)
    }

    fn total_sum(&self) -> u32 {
        self.prefix_sum(self.size)
    }
}

/// Genomic variant record
#[derive(Debug, Clone)]
struct GenotypeRecord {
    chromosome: String,
    position: usize,
    genotype_count: u8,
    allele_freq: Option<f64>,
}

/// Statistics collector for genomic data
#[derive(Debug, Default)]
struct GenomicStats {
    total_variants: usize,
    chromosomes: HashMap<String, usize>,
    min_position: usize,
    max_position: usize,
    total_alleles: u64,
}

impl GenomicStats {
    fn update(&mut self, record: &GenotypeRecord) {
        self.total_variants += 1;
        *self.chromosomes.entry(record.chromosome.clone()).or_insert(0) += 1;
        
        if self.min_position == 0 || record.position < self.min_position {
            self.min_position = record.position;
        }
        if record.position > self.max_position {
            self.max_position = record.position;
        }
        
        self.total_alleles += record.genotype_count as u64;
    }

    fn print_summary(&self) {
        info!("=== Genomic Data Summary ===");
        info!("Total variants: {}", self.total_variants);
        info!("Total alleles: {}", self.total_alleles);
        info!("Position range: {} - {}", self.min_position, self.max_position);
        info!("Chromosomes found:");
        for (chr, count) in &self.chromosomes {
            info!("  {}: {} variants", chr, count);
        }
    }
}

/// Parse a single line from the genotype file
fn parse_genotype_line(line: &str, line_num: usize) -> Result<GenotypeRecord> {
    let tokens: Vec<&str> = line.trim().split_whitespace().collect();
    
    if tokens.len() < 3 {
        return Err(anyhow::anyhow!(
            "Line {}: Expected at least 3 columns, found {}. Line: '{}'", 
            line_num, tokens.len(), line
        ));
    }

    let chromosome = tokens[0].to_string();
    let position = tokens[1].parse::<usize>()
        .with_context(|| format!("Line {}: Invalid position '{}'", line_num, tokens[1]))?;
    let genotype_count = tokens[2].parse::<u8>()
        .with_context(|| format!("Line {}: Invalid genotype count '{}'", line_num, tokens[2]))?;

    // Optional allele frequency in 4th column
    let allele_freq = if tokens.len() > 3 {
        Some(tokens[3].parse::<f64>()
            .with_context(|| format!("Line {}: Invalid allele frequency '{}'", line_num, tokens[3]))?)
    } else {
        None
    };

    Ok(GenotypeRecord {
        chromosome,
        position,
        genotype_count,
        allele_freq,
    })
}

/// Load genomic data from file with progress tracking
async fn load_genomic_data(file_path: &str) -> Result<(Vec<GenotypeRecord>, GenomicStats)> {
    info!("Loading genomic data from: {}", file_path);
    
    let file = File::open(file_path)
        .with_context(|| format!("Failed to open file: {}", file_path))?;
    
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut stats = GenomicStats::default();
    
    // Count lines for progress bar
    let total_lines = BufReader::new(File::open(file_path)?)
        .lines()
        .count();
    
    let pb = ProgressBar::new(total_lines as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"));

    let mut line_num = 0;
    let mut error_count = 0;
    
    for line_result in reader.lines() {
        line_num += 1;
        pb.inc(1);
        
        let line = line_result
            .with_context(|| format!("Failed to read line {}", line_num))?;
        
        // Skip empty lines and comments
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        
        match parse_genotype_line(&line, line_num) {
            Ok(record) => {
                stats.update(&record);
                records.push(record);
            }
            Err(e) => {
                error_count += 1;
                if error_count <= 10 {
                    warn!("Parse error: {}", e);
                } else if error_count == 11 {
                    warn!("Suppressing further parse errors...");
                }
            }
        }
    }
    
    pb.finish_with_message("Data loading complete");
    
    if error_count > 0 {
        warn!("Total parse errors: {}", error_count);
    }
    
    stats.print_summary();
    Ok((records, stats))
}

/// Build Fenwick tree from genomic records
fn build_fenwick_tree(records: &[GenotypeRecord], max_index: usize) -> Result<FenwickTree> {
    info!("Building Fenwick tree with max index: {}", max_index);
    
    let mut tree = FenwickTree::new(max_index);
    let pb = ProgressBar::new(records.len() as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} Building index...")
        .unwrap());

    for record in records {
        if record.position <= max_index {
            tree.update(record.position, record.genotype_count as u32)
                .with_context(|| format!("Failed to update tree at position {}", record.position))?;
        } else {
            warn!("Position {} exceeds max index {}, skipping", record.position, max_index);
        }
        pb.inc(1);
    }
    
    pb.finish_with_message("Fenwick tree built successfully");
    
    info!("Total alleles indexed: {}", tree.total_sum());
    Ok(tree)
}

/// Perform parallel range queries
async fn perform_queries(tree: Arc<FenwickTree>, num_tasks: usize, max_index: usize) -> Result<Vec<(usize, usize, u32)>> {
    info!("Performing {} parallel range queries", num_tasks);
    
    let mut handles = Vec::new();
    let results = Arc::new(Mutex::new(Vec::new()));
    
    for i in 0..num_tasks {
        let tree_clone = Arc::clone(&tree);
        let results_clone = Arc::clone(&results);
        
        handles.push(task::spawn(async move {
            let chunk_size = max_index / num_tasks;
            let start = i * chunk_size + 1;
            let end = if i == num_tasks - 1 { max_index } else { (i + 1) * chunk_size };
            
            match tree_clone.range_sum(start, end) {
                Ok(sum_val) => {
                    let mut results_guard = results_clone.lock().unwrap();
                    results_guard.push((start, end, sum_val));
                    info!("Task {}: Range sum from {} to {} is {}", i + 1, start, end, sum_val);
                }
                Err(e) => {
                    error!("Task {}: Query failed: {}", i + 1, e);
                }
            }
        }));
    }

    // Wait for all tasks to complete
    for handle in handles {
        handle.await
            .with_context(|| "Failed to complete query task")?;
    }
    
    let results = Arc::try_unwrap(results).unwrap().into_inner().unwrap();
    Ok(results)
}

/// Write results to output file
fn write_results(
    output_path: &str, 
    query_results: &[(usize, usize, u32)], 
    stats: &GenomicStats,
    tree: &FenwickTree
) -> Result<()> {
    info!("Writing results to: {}", output_path);
    
    let mut file = File::create(output_path)
        .with_context(|| format!("Failed to create output file: {}", output_path))?;
    
    writeln!(file, "# Genomic Sparse Index Results")?;
    writeln!(file, "# Generated by rust_sparse_index")?;
    writeln!(file, "# Total variants: {}", stats.total_variants)?;
    writeln!(file, "# Total alleles: {}", stats.total_alleles)?;
    writeln!(file, "# Position range: {} - {}", stats.min_position, stats.max_position)?;
    writeln!(file, "# Tree total sum: {}", tree.total_sum())?;
    writeln!(file, "#")?;
    writeln!(file, "# Range Query Results:")?;
    writeln!(file, "start_pos\tend_pos\tallele_sum")?;
    
    for (start, end, sum) in query_results {
        writeln!(file, "{}\t{}\t{}", start, end, sum)?;
    }
    
    Ok(())
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = Args::parse();
    
    // Initialize logger
    if args.verbose {
        env_logger::Builder::from_default_env()
            .filter_level(log::LevelFilter::Debug)
            .init();
    } else {
        env_logger::Builder::from_default_env()
            .filter_level(log::LevelFilter::Info)
            .init();
    }
    
    info!("Starting genomic sparse index builder");
    info!("Input file: {}", args.input);
    info!("Output file: {}", args.output);
    info!("Max index: {}", args.max_index);
    info!("Number of query tasks: {}", args.num_tasks);
    
    // Validate input file exists
    if !Path::new(&args.input).exists() {
        return Err(anyhow::anyhow!("Input file does not exist: {}", args.input));
    }
    
    // Load genomic data
    let (records, stats) = load_genomic_data(&args.input).await?;
    
    if records.is_empty() {
        return Err(anyhow::anyhow!("No valid records found in input file"));
    }
    
    // Build Fenwick tree
    let tree = build_fenwick_tree(&records, args.max_index)?;
    let tree_arc = Arc::new(tree);
    
    // Perform parallel queries
    let query_results = perform_queries(Arc::clone(&tree_arc), args.num_tasks, args.max_index).await?;
    
    // Write results
    write_results(&args.output, &query_results, &stats, &tree_arc)?;
    
    info!("Processing complete! Results written to: {}", args.output);
    
    Ok(())
}