use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::collections::HashMap;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn, error, debug};
use rayon::prelude::*;
use ndarray::{Array2, Array1, s, Axis};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use rand_distr::{Distribution, Normal};
use serde::{Deserialize, Serialize};
use statrs::statistics::{Statistics, Data};

/// Genomic Algorithms Tool for Allele Frequency, LD Analysis, and HMM Phasing
#[derive(Parser, Debug)]
#[command(name = "rust_core_algorithms")]
#[command(about = "Advanced genomic algorithms for population genetics analysis")]
struct Args {
    /// Input genotype CSV file
    #[arg(short, long)]
    input: String,
    
    /// Output results file
    #[arg(short, long, default_value = "results.txt")]
    output: String,
    
    /// Output format (txt, json, csv)
    #[arg(short, long, default_value = "txt")]
    format: String,
    
    /// Number of HMM states for phasing
    #[arg(long, default_value = "3")]
    hmm_states: usize,
    
    /// Random seed for reproducibility
    #[arg(long, default_value = "42")]
    seed: u64,
    
    /// Enable LD calculation (can be slow for large datasets)
    #[arg(long)]
    compute_ld: bool,
    
    /// Maximum number of variant pairs for LD computation
    #[arg(long, default_value = "1000")]
    max_ld_pairs: usize,
    
    /// Enable HMM phasing for all individuals
    #[arg(long)]
    phase_all: bool,
    
    /// Minimum allele frequency threshold
    #[arg(long, default_value = "0.01")]
    min_maf: f64,
    
    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,
}

/// Genomic variant information
#[derive(Debug, Clone, Serialize, Deserialize)]
struct Variant {
    id: String,
    chromosome: String,
    position: u64,
    ref_allele: String,
    alt_allele: String,
    allele_frequency: f64,
    maf: f64,
    call_rate: f64,
}

/// Linkage disequilibrium result
#[derive(Debug, Clone, Serialize, Deserialize)]
struct LDResult {
    variant1: String,
    variant2: String,
    r_squared: f64,
    d_prime: f64,
    distance: i64,
}

/// HMM phasing result
#[derive(Debug, Clone, Serialize, Deserialize)]
struct PhasingResult {
    individual_id: String,
    haplotype1: Vec<u8>,
    haplotype2: Vec<u8>,
    confidence_scores: Vec<f64>,
    switch_points: Vec<usize>,
}

/// Comprehensive genomic analysis results
#[derive(Debug, Serialize, Deserialize)]
struct GenomicResults {
    summary: AnalysisSummary,
    variants: Vec<Variant>,
    ld_results: Vec<LDResult>,
    phasing_results: Vec<PhasingResult>,
    processing_time_ms: u128,
}

#[derive(Debug, Serialize, Deserialize)]
struct AnalysisSummary {
    total_variants: usize,
    total_individuals: usize,
    mean_call_rate: f64,
    mean_maf: f64,
    variants_passing_qc: usize,
    ld_pairs_computed: usize,
    individuals_phased: usize,
}

/// Enhanced Hidden Markov Model for haplotype phasing
#[derive(Debug, Clone)]
struct GenomicHMM {
    transition_matrix: Array2<f64>,
    emission_matrix: Array2<f64>,
    initial_probs: Array1<f64>,
    num_states: usize,
    num_emissions: usize,
}

impl GenomicHMM {
    fn new(num_states: usize, num_emissions: usize, seed: u64) -> Result<Self> {
        let mut rng = StdRng::seed_from_u64(seed);
        
        // Initialize with biologically reasonable parameters
        let mut transition_matrix = Array2::zeros((num_states, num_states));
        let mut emission_matrix = Array2::zeros((num_states, num_emissions));
        let mut initial_probs = Array1::zeros(num_states);
        
        // Set transition probabilities (higher probability for staying in same state)
        for i in 0..num_states {
            for j in 0..num_states {
                if i == j {
                    transition_matrix[[i, j]] = 0.95 + rng.gen::<f64>() * 0.04; // High self-transition
                } else {
                    transition_matrix[[i, j]] = 0.001 + rng.gen::<f64>() * 0.009; // Low switch probability
                }
            }
        }
        
        // Set emission probabilities based on genotype patterns
        for i in 0..num_states {
            for j in 0..num_emissions {
                emission_matrix[[i, j]] = match (i, j) {
                    (0, 0) => 0.8 + rng.gen::<f64>() * 0.15, // Homozygous ref state emits mostly 0
                    (0, 1) => 0.15 + rng.gen::<f64>() * 0.05,
                    (0, 2) => 0.05 + rng.gen::<f64>() * 0.05,
                    (1, 0) => 0.1 + rng.gen::<f64>() * 0.1,   // Heterozygous state
                    (1, 1) => 0.8 + rng.gen::<f64>() * 0.15,
                    (1, 2) => 0.1 + rng.gen::<f64>() * 0.1,
                    (2, 0) => 0.05 + rng.gen::<f64>() * 0.05, // Homozygous alt state emits mostly 2
                    (2, 1) => 0.15 + rng.gen::<f64>() * 0.05,
                    (2, 2) => 0.8 + rng.gen::<f64>() * 0.15,
                    _ => 0.33 + rng.gen::<f64>() * 0.34,      // Default equal probability
                };
            }
        }
        
        // Initialize state probabilities
        for i in 0..num_states {
            initial_probs[i] = 1.0 / num_states as f64 + rng.gen::<f64>() * 0.1;
        }
        
        // Normalize all probability matrices
        Self::normalize_matrix(&mut transition_matrix, Axis(1))?;
        Self::normalize_matrix(&mut emission_matrix, Axis(1))?;
        initial_probs /= initial_probs.sum();
        
        Ok(GenomicHMM {
            transition_matrix,
            emission_matrix,
            initial_probs,
            num_states,
            num_emissions,
        })
    }
    
    fn normalize_matrix(matrix: &mut Array2<f64>, axis: Axis) -> Result<()> {
        for mut row in matrix.lanes_mut(axis) {
            let sum = row.sum();
            if sum > 0.0 {
                row /= sum;
            } else {
                // Handle zero-sum rows by setting uniform distribution
                row.fill(1.0 / row.len() as f64);
            }
        }
        Ok(())
    }
    
    fn viterbi_decode(&self, observations: &[usize]) -> Result<(Vec<usize>, f64)> {
        if observations.is_empty() {
            return Ok((Vec::new(), 0.0));
        }
        
        let n_obs = observations.len();
        let mut dp = Array2::<f64>::zeros((n_obs, self.num_states));
        let mut backtrack = Array2::<usize>::zeros((n_obs, self.num_states));
        
        // Initialization step
        for s in 0..self.num_states {
            if observations[0] < self.num_emissions {
                let emission_prob = self.emission_matrix[[s, observations[0]]];
                dp[[0, s]] = (self.initial_probs[s] * emission_prob).ln();
            } else {
                return Err(anyhow::anyhow!("Invalid observation: {}", observations[0]));
            }
        }
        
        // Recursion step (using log probabilities for numerical stability)
        for t in 1..n_obs {
            if observations[t] >= self.num_emissions {
                return Err(anyhow::anyhow!("Invalid observation at position {}: {}", t, observations[t]));
            }
            
            for s in 0..self.num_states {
                let emission_prob = self.emission_matrix[[s, observations[t]]];
                let emission_log = emission_prob.ln();
                
                let (max_prev_idx, max_val) = (0..self.num_states)
                    .map(|s_prev| {
                        let transition_prob = self.transition_matrix[[s_prev, s]];
                        let score = dp[[t-1, s_prev]] + transition_prob.ln() + emission_log;
                        (s_prev, score)
                    })
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap();
                
                dp[[t, s]] = max_val;
                backtrack[[t, s]] = max_prev_idx;
            }
        }
        
        // Termination step
        let (max_end_state, max_prob) = (0..self.num_states)
            .map(|s| (s, dp[[n_obs-1, s]]))
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap();
        
        // Path reconstruction
        let mut path = vec![0; n_obs];
        path[n_obs-1] = max_end_state;
        let mut current_state = max_end_state;
        
        for t in (0..n_obs-1).rev() {
            current_state = backtrack[[t+1, current_state]];
            path[t] = current_state;
        }
        
        Ok((path, max_prob))
    }
    
    fn forward_backward(&self, observations: &[usize]) -> Result<(Array2<f64>, f64)> {
        let n_obs = observations.len();
        let mut forward = Array2::<f64>::zeros((n_obs, self.num_states));
        let mut backward = Array2::<f64>::zeros((n_obs, self.num_states));
        
        // Forward algorithm
        for s in 0..self.num_states {
            forward[[0, s]] = self.initial_probs[s] * self.emission_matrix[[s, observations[0]]];
        }
        
        for t in 1..n_obs {
            for s in 0..self.num_states {
                let mut sum = 0.0;
                for s_prev in 0..self.num_states {
                    sum += forward[[t-1, s_prev]] * self.transition_matrix[[s_prev, s]];
                }
                forward[[t, s]] = sum * self.emission_matrix[[s, observations[t]]];
            }
        }
        
        // Backward algorithm
        for s in 0..self.num_states {
            backward[[n_obs-1, s]] = 1.0;
        }
        
        for t in (0..n_obs-1).rev() {
            for s in 0..self.num_states {
                let mut sum = 0.0;
                for s_next in 0..self.num_states {
                    sum += self.transition_matrix[[s, s_next]] 
                         * self.emission_matrix[[s_next, observations[t+1]]] 
                         * backward[[t+1, s_next]];
                }
                backward[[t, s]] = sum;
            }
        }
        
        // Calculate total likelihood
        let likelihood = forward.slice(s![n_obs-1, ..]).sum();
        
        // Calculate posterior probabilities
        let mut posterior = Array2::<f64>::zeros((n_obs, self.num_states));
        for t in 0..n_obs {
            for s in 0..self.num_states {
                posterior[[t, s]] = forward[[t, s]] * backward[[t, s]] / likelihood;
            }
        }
        
        Ok((posterior, likelihood))
    }
}

/// Compute allele frequencies with quality metrics
fn compute_allele_frequencies(genotypes: &Array2<u8>) -> Result<Vec<Variant>> {
    let (n_individuals, n_variants) = genotypes.dim();
    info!("Computing allele frequencies for {} variants across {} individuals", n_variants, n_individuals);
    
    let pb = ProgressBar::new(n_variants as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} Computing frequencies...")
        .unwrap());
    
    let variants: Vec<Variant> = (0..n_variants).into_par_iter().map(|var_idx| {
        pb.inc(1);
        
        let genotype_column = genotypes.column(var_idx);
        let mut allele_counts = [0u32; 3]; // 0, 1, 2 copies of alt allele
        let mut missing_count = 0u32;
        
        for &genotype in genotype_column.iter() {
            match genotype {
                0 | 1 | 2 => allele_counts[genotype as usize] += 1,
                _ => missing_count += 1,
            }
        }
        
        let total_called = allele_counts.iter().sum::<u32>();
        let call_rate = total_called as f64 / n_individuals as f64;
        
        let total_alleles = (allele_counts[0] * 0 + allele_counts[1] * 1 + allele_counts[2] * 2) as f64;
        let max_alleles = total_called as f64 * 2.0;
        
        let allele_frequency = if max_alleles > 0.0 {
            total_alleles / max_alleles
        } else {
            0.0
        };
        
        let maf = allele_frequency.min(1.0 - allele_frequency);
        
        Variant {
            id: format!("variant_{}", var_idx + 1),
            chromosome: format!("chr{}", ((var_idx / 10000) % 22) + 1),
            position: (var_idx * 1000 + 1000) as u64,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            allele_frequency,
            maf,
            call_rate,
        }
    }).collect();
    
    pb.finish_with_message("Allele frequency computation completed");
    Ok(variants)
}

/// Compute linkage disequilibrium between variant pairs
fn compute_linkage_disequilibrium(
    genotypes: &Array2<u8>, 
    variants: &[Variant], 
    max_pairs: usize
) -> Result<Vec<LDResult>> {
    let n_variants = variants.len();
    info!("Computing LD for up to {} variant pairs", max_pairs);
    
    let total_pairs = std::cmp::min(max_pairs, (n_variants * (n_variants - 1)) / 2);
    let pb = ProgressBar::new(total_pairs as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} Computing LD...")
        .unwrap());
    
    let mut ld_results = Vec::new();
    let mut pairs_computed = 0;
    
    'outer: for i in 0..n_variants {
        for j in (i+1)..n_variants {
            if pairs_computed >= max_pairs {
                break 'outer;
            }
            
            pb.inc(1);
            pairs_computed += 1;
            
            let geno_i = genotypes.column(i);
            let geno_j = genotypes.column(j);
            
            // Count haplotype frequencies
            let mut counts = [[0u32; 3]; 3];
            let mut total_count = 0u32;
            
            for (&gi, &gj) in geno_i.iter().zip(geno_j.iter()) {
                if gi <= 2 && gj <= 2 {
                    counts[gi as usize][gj as usize] += 1;
                    total_count += 1;
                }
            }
            
            if total_count < 10 {
                continue; // Skip pairs with insufficient data
            }
            
            // Calculate allele frequencies
            let p_a = (counts[1][0] + counts[1][1] + counts[1][2] + 2 * counts[2][0] + 2 * counts[2][1] + 2 * counts[2][2]) as f64 / (2.0 * total_count as f64);
            let p_b = (counts[0][1] + counts[1][1] + counts[2][1] + 2 * counts[0][2] + 2 * counts[1][2] + 2 * counts[2][2]) as f64 / (2.0 * total_count as f64);
            
            // Calculate haplotype frequency for AB
            let p_ab = (counts[1][1] + 2.0 * counts[2][2] as f64 + counts[1][2] as f64 + counts[2][1] as f64) / (2.0 * total_count as f64);
            
            // Calculate D and r^2
            let d = p_ab - (p_a * p_b);
            let d_max = if d > 0.0 {
                ((1.0 - p_a) * (1.0 - p_b)).min(p_a * p_b)
            } else {
                (p_a * (1.0 - p_b)).min((1.0 - p_a) * p_b)
            };
            
            let d_prime = if d_max != 0.0 { d.abs() / d_max } else { 0.0 };
            
            let r_squared = if p_a > 0.0 && p_a < 1.0 && p_b > 0.0 && p_b < 1.0 {
                (d * d) / (p_a * (1.0 - p_a) * p_b * (1.0 - p_b))
            } else {
                0.0
            };
            
            let distance = (variants[j].position as i64) - (variants[i].position as i64);
            
            ld_results.push(LDResult {
                variant1: variants[i].id.clone(),
                variant2: variants[j].id.clone(),
                r_squared: r_squared.max(0.0).min(1.0),
                d_prime: d_prime.max(0.0).min(1.0),
                distance,
            });
        }
    }
    
    pb.finish_with_message("LD computation completed");
    Ok(ld_results)
}

/// Perform HMM-based phasing for individuals
fn perform_hmm_phasing(
    genotypes: &Array2<u8>, 
    hmm_states: usize,
    seed: u64,
    phase_all: bool
) -> Result<Vec<PhasingResult>> {
    let (n_individuals, n_variants) = genotypes.dim();
    let individuals_to_phase = if phase_all { n_individuals } else { std::cmp::min(5, n_individuals) };
    
    info!("Performing HMM phasing for {} individuals with {} states", individuals_to_phase, hmm_states);
    
    let pb = ProgressBar::new(individuals_to_phase as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} Phasing individuals...")
        .unwrap());
    
    let phasing_results: Result<Vec<PhasingResult>> = (0..individuals_to_phase).into_par_iter().map(|ind_idx| {
        pb.inc(1);
        
        let individual_genotypes = genotypes.row(ind_idx);
        let observations: Vec<usize> = individual_genotypes.iter()
            .map(|&g| std::cmp::min(g as usize, 2))
            .collect();
        
        let hmm = GenomicHMM::new(hmm_states, 3, seed + ind_idx as u64)?;
        let (path, likelihood) = hmm.viterbi_decode(&observations)?;
        let (posterior, _) = hmm.forward_backward(&observations)?;
        
        // Extract confidence scores and identify switch points
        let confidence_scores: Vec<f64> = (0..n_variants).map(|t| {
            posterior.row(t).iter().fold(0.0, |max_val, &val| max_val.max(val))
        }).collect();
        
        let mut switch_points = Vec::new();
        for t in 1..path.len() {
            if path[t] != path[t-1] && confidence_scores[t] > 0.8 {
                switch_points.push(t);
            }
        }
        
        // Generate haplotypes based on HMM path and genotypes
        let mut haplotype1 = Vec::new();
        let mut haplotype2 = Vec::new();
        
        for (t, &genotype) in observations.iter().enumerate() {
            match (genotype, path[t]) {
                (0, _) => { haplotype1.push(0); haplotype2.push(0); },
                (2, _) => { haplotype1.push(1); haplotype2.push(1); },
                (1, state) => {
                    if state % 2 == 0 {
                        haplotype1.push(0); haplotype2.push(1);
                    } else {
                        haplotype1.push(1); haplotype2.push(0);
                    }
                },
                _ => { haplotype1.push(0); haplotype2.push(0); }, // Default for missing data
            }
        }
        
        Ok(PhasingResult {
            individual_id: format!("individual_{}", ind_idx + 1),
            haplotype1,
            haplotype2,
            confidence_scores,
            switch_points,
        })
    }).collect();
    
    pb.finish_with_message("HMM phasing completed");
    phasing_results
}

/// Load genotype data from CSV file
fn load_genotype_data(file_path: &str) -> Result<Array2<u8>> {
    info!("Loading genotype data from: {}", file_path);
    
    let file = File::open(file_path)
        .with_context(|| format!("Failed to open file: {}", file_path))?;
    
    let reader = BufReader::new(file);
    let mut genotypes: Vec<Vec<u8>> = Vec::new();
    let mut line_count = 0;
    
    for line_result in reader.lines() {
        line_count += 1;
        let line = line_result
            .with_context(|| format!("Failed to read line {}", line_count))?;
        
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        
        let row: Result<Vec<u8>> = line.split(',')
            .enumerate()
            .map(|(col_idx, val)| {
                val.trim().parse::<u8>()
                    .with_context(|| format!("Invalid genotype at line {}, column {}: '{}'", line_count, col_idx + 1, val))
            })
            .collect();
        
        match row {
            Ok(genotype_row) => genotypes.push(genotype_row),
            Err(e) => {
                warn!("Skipping line {}: {}", line_count, e);
                continue;
            }
        }
    }
    
    if genotypes.is_empty() {
        return Err(anyhow::anyhow!("No valid genotype data found in file"));
    }
    
    let n_individuals = genotypes.len();
    let n_variants = genotypes[0].len();
    
    // Validate that all rows have the same length
    for (i, row) in genotypes.iter().enumerate() {
        if row.len() != n_variants {
            return Err(anyhow::anyhow!(
                "Row {} has {} variants, expected {}",
                i + 1, row.len(), n_variants
            ));
        }
    }
    
    // Convert to ndarray
    let flat_data: Vec<u8> = genotypes.into_iter().flatten().collect();
    let genotype_matrix = Array2::from_shape_vec((n_individuals, n_variants), flat_data)
        .context("Failed to create genotype matrix")?;
    
    info!("Loaded {} individuals × {} variants", n_individuals, n_variants);
    Ok(genotype_matrix)
}

/// Write results to file in specified format
fn write_results(results: &GenomicResults, output_path: &str, format: &str) -> Result<()> {
    info!("Writing results to: {} (format: {})", output_path, format);
    
    match format.to_lowercase().as_str() {
        "json" => {
            let json_output = serde_json::to_string_pretty(results)
                .context("Failed to serialize results to JSON")?;
            std::fs::write(output_path, json_output)
                .context("Failed to write JSON results")?;
        },
        "csv" => {
            let mut file = File::create(output_path)
                .context("Failed to create CSV output file")?;
            
            // Write variants CSV format
            writeln!(file, "variant_id,chromosome,position,allele_frequency,maf,call_rate")?;
            for variant in &results.variants {
                writeln!(file, "{},{},{},{:.6},{:.6},{:.4}", 
                    variant.id, variant.chromosome, variant.position,
                    variant.allele_frequency, variant.maf, variant.call_rate)?;
            }
        },
        _ => {
            // Default text format
            let mut file = File::create(output_path)
                .context("Failed to create text output file")?;
            
            writeln!(file, "# Genomic Analysis Results")?;
            writeln!(file, "# Processing time: {} ms", results.processing_time_ms)?;
            writeln!(file, "#")?;
            
            writeln!(file, "## Summary Statistics")?;
            writeln!(file, "Total variants: {}", results.summary.total_variants)?;
            writeln!(file, "Total individuals: {}", results.summary.total_individuals)?;
            writeln!(file, "Mean call rate: {:.4}", results.summary.mean_call_rate)?;
            writeln!(file, "Mean MAF: {:.4}", results.summary.mean_maf)?;
            writeln!(file, "Variants passing QC: {}", results.summary.variants_passing_qc)?;
            writeln!(file, "LD pairs computed: {}", results.summary.ld_pairs_computed)?;
            writeln!(file, "Individuals phased: {}", results.summary.individuals_phased)?;
            writeln!(file)?;
            
            writeln!(file, "## Variant Information (first 20)")?;
            writeln!(file, "ID\tChromosome\tPosition\tAF\tMAF\tCall_Rate")?;
            for variant in results.variants.iter().take(20) {
                writeln!(file, "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.4}",
                    variant.id, variant.chromosome, variant.position,
                    variant.allele_frequency, variant.maf, variant.call_rate)?;
            }
            
            if !results.ld_results.is_empty() {
                writeln!(file)?;
                writeln!(file, "## Linkage Disequilibrium Results (first 20)")?;
                writeln!(file, "Variant1\tVariant2\tR_squared\tD_prime\tDistance")?;
                for ld in results.ld_results.iter().take(20) {
                    writeln!(file, "{}\t{}\t{:.6}\t{:.6}\t{}",
                        ld.variant1, ld.variant2, ld.r_squared, ld.d_prime, ld.distance)?;
                }
            }
            
            if !results.phasing_results.is_empty() {
                writeln!(file)?;
                writeln!(file, "## Phasing Results Summary")?;
                for phasing in results.phasing_results.iter().take(5) {
                    writeln!(file, "Individual: {}", phasing.individual_id)?;
                    writeln!(file, "  Switch points: {:?}", phasing.switch_points)?;
                    writeln!(file, "  Mean confidence: {:.4}", 
                        phasing.confidence_scores.iter().sum::<f64>() / phasing.confidence_scores.len() as f64)?;
                    writeln!(file, "  Haplotype 1 (first 20): {:?}", 
                        &phasing.haplotype1.iter().take(20).collect::<Vec<_>>())?;
                    writeln!(file, "  Haplotype 2 (first 20): {:?}", 
                        &phasing.haplotype2.iter().take(20).collect::<Vec<_>>())?;
                    writeln!(file)?;
                }
            }
        }
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
    
    let start_time = Instant::now();
    
    info!("Starting genomic analysis pipeline");
    info!("Input file: {}", args.input);
    info!("Output file: {}", args.output);
    info!("Output format: {}", args.format);
    info!("HMM states: {}", args.hmm_states);
    info!("Random seed: {}", args.seed);
    info!("Compute LD: {}", args.compute_ld);
    info!("Phase all individuals: {}", args.phase_all);
    info!("Minimum MAF: {}", args.min_maf);
    
    // Validate input file
    if !Path::new(&args.input).exists() {
        return Err(anyhow::anyhow!("Input file does not exist: {}", args.input));
    }
    
    // Load genotype data
    let genotypes = load_genotype_data(&args.input)?;
    let (n_individuals, n_variants) = genotypes.dim();
    
    if n_variants == 0 || n_individuals == 0 {
        return Err(anyhow::anyhow!("No valid genotype data found"));
    }
    
    // Compute allele frequencies
    let variants = compute_allele_frequencies(&genotypes)?;
    
    // Filter variants by MAF
    let variants_passing_qc = variants.iter()
        .filter(|v| v.maf >= args.min_maf && v.call_rate >= 0.95)
        .count();
    
    info!("{} variants pass QC filters (MAF >= {}, call rate >= 0.95)", 
          variants_passing_qc, args.min_maf);
    
    // Compute linkage disequilibrium if requested
    let ld_results = if args.compute_ld {
        compute_linkage_disequilibrium(&genotypes, &variants, args.max_ld_pairs)?
    } else {
        info!("Skipping LD computation (use --compute-ld to enable)");
        Vec::new()
    };
    
    // Perform HMM phasing
    let phasing_results = perform_hmm_phasing(
        &genotypes, 
        args.hmm_states, 
        args.seed, 
        args.phase_all
    )?;
    
    // Calculate summary statistics
    let mean_call_rate = variants.iter().map(|v| v.call_rate).sum::<f64>() / variants.len() as f64;
    let mean_maf = variants.iter().map(|v| v.maf).sum::<f64>() / variants.len() as f64;
    
    let processing_time = start_time.elapsed();
    
    let results = GenomicResults {
        summary: AnalysisSummary {
            total_variants: n_variants,
            total_individuals: n_individuals,
            mean_call_rate,
            mean_maf,
            variants_passing_qc,
            ld_pairs_computed: ld_results.len(),
            individuals_phased: phasing_results.len(),
        },
        variants,
        ld_results,
        phasing_results,
        processing_time_ms: processing_time.as_millis(),
    };
    
    // Write results
    write_results(&results, &args.output, &args.format)?;
    
    info!("Analysis completed successfully!");
    info!("Processing time: {:.2} seconds", processing_time.as_secs_f64());
    info!("Results written to: {}", args.output);
    
    // Print summary to console
    println!("\n=== Analysis Summary ===");
    println!("Variants processed: {}", results.summary.total_variants);
    println!("Individuals analyzed: {}", results.summary.total_individuals);
    println!("Mean call rate: {:.4}", results.summary.mean_call_rate);
    println!("Mean MAF: {:.4}", results.summary.mean_maf);
    println!("Variants passing QC: {}", results.summary.variants_passing_qc);
    println!("LD pairs computed: {}", results.summary.ld_pairs_computed);
    println!("Individuals phased: {}", results.summary.individuals_phased);
    println!("Processing time: {:.2} seconds", processing_time.as_secs_f64());
    
    Ok(())
}
            