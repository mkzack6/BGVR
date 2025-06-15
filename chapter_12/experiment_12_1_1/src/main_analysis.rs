use csv::{Reader, Writer};
use ndarray::{Array1, Array2, Axis};
use ndarray_linalg::Solve;
use ndarray_stats::QuantileExt;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use clap::Parser;
use anyhow::{Context, Result};
use thiserror::Error;
use rayon::prelude::*;
use log::{info, warn, error};

#[derive(Error, Debug)]
pub enum AnalysisError {
    #[error("Data parsing error: {0}")]
    DataParsing(String),
    #[error("Matrix operation error: {0}")]
    MatrixOperation(String),
    #[error("File I/O error: {0}")]
    FileIO(String),
}

/// Population genomics analysis pipeline
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input genotype file
    #[arg(short, long, default_value = "data/genotypes.csv")]
    genotypes: String,
    
    /// Sample metadata file
    #[arg(short, long, default_value = "data/sample_metadata.csv")]
    metadata: String,
    
    /// Output directory
    #[arg(short, long, default_value = "results")]
    output: String,
    
    /// Minor allele frequency threshold
    #[arg(long, default_value_t = 0.05)]
    maf_threshold: f64,
    
    /// Missing data threshold per variant
    #[arg(long, default_value_t = 0.1)]
    missing_threshold: f64,
    
    /// Missing data threshold per sample
    #[arg(long, default_value_t = 0.1)]
    sample_missing_threshold: f64,
    
    /// Number of principal components to compute
    #[arg(long, default_value_t = 10)]
    num_pcs: usize,
    
    /// Perform Hardy-Weinberg equilibrium test
    #[arg(long)]
    hwe_test: bool,
    
    /// Perform population structure analysis
    #[arg(long)]
    population_structure: bool,
}

#[derive(Debug, Clone)]
struct SampleInfo {
    id: String,
    population: String,
    sex: String,
    age: i32,
    phenotype: f64,
}

#[derive(Debug)]
struct QualityMetrics {
    total_variants: usize,
    variants_after_maf: usize,
    variants_after_missing: usize,
    total_samples: usize,
    samples_after_qc: usize,
    mean_missing_per_sample: f64,
    mean_missing_per_variant: f64,
}

#[tokio::main]
async fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();
    
    info!("Starting population genomics analysis pipeline");
    info!("Parameters: MAF threshold = {:.3}, Missing threshold = {:.3}", 
          args.maf_threshold, args.missing_threshold);
    
    // Create output directory
    std::fs::create_dir_all(&args.output)?;
    
    // Load sample metadata
    let sample_info = load_sample_metadata(&args.metadata)
        .context("Failed to load sample metadata")?;
    info!("Loaded metadata for {} samples", sample_info.len());
    
    // Load and process genotype data
    let (genotype_matrix, sample_ids, variant_ids, qc_metrics) = 
        load_and_filter_genotypes(&args, &sample_info)
        .context("Failed to load and filter genotype data")?;
    
    // Save QC report
    save_qc_report(&args.output, &qc_metrics)?;
    
    // Perform PCA
    info!("Performing Principal Component Analysis...");
    let pca_results = perform_pca(&genotype_matrix, args.num_pcs)
        .context("Failed to perform PCA")?;
    
    // Save PCA results
    save_pca_results(&args.output, &pca_results, &sample_ids, &sample_info)?;
    
    // Calculate allele frequencies per population
    if args.population_structure {
        info!("Calculating population-specific allele frequencies...");
        calculate_population_structure(&genotype_matrix, &sample_ids, &sample_info, &args.output)?;
    }
    
    // Hardy-Weinberg equilibrium testing
    if args.hwe_test {
        info!("Performing Hardy-Weinberg equilibrium tests...");
        perform_hwe_tests(&genotype_matrix, &variant_ids, &args.output)?;
    }
    
    info!("Analysis complete! Results saved to {}", args.output);
    Ok(())
}

fn load_sample_metadata(metadata_path: &str) -> Result<HashMap<String, SampleInfo>> {
    let mut sample_info = HashMap::new();
    let mut rdr = Reader::from_path(metadata_path)
        .context("Failed to open metadata file")?;
    
    for result in rdr.records() {
        let record = result.context("Failed to read metadata record")?;
        
        if record.len() < 5 {
            warn!("Skipping incomplete metadata record: {:?}", record);
            continue;
        }
        
        let info = SampleInfo {
            id: record[0].to_string(),
            population: record[1].to_string(),
            sex: record[2].to_string(),
            age: record[3].parse().unwrap_or(0),
            phenotype: record[4].parse().unwrap_or(0.0),
        };
        
        sample_info.insert(info.id.clone(), info);
    }
    
    Ok(sample_info)
}

fn load_and_filter_genotypes(
    args: &Args, 
    sample_info: &HashMap<String, SampleInfo>
) -> Result<(Array2<f64>, Vec<String>, Vec<String>, QualityMetrics)> {
    
    info!("Loading genotype data from {}", args.genotypes);
    let mut rdr = Reader::from_path(&args.genotypes)
        .context("Failed to open genotype file")?;
    
    // Read header to get variant IDs
    let headers = rdr.headers()
        .context("Failed to read genotype file headers")?;
    let variant_ids: Vec<String> = headers.iter().skip(1).map(|s| s.to_string()).collect();
    let total_variants = variant_ids.len();
    
    info!("Found {} variants", total_variants);
    
    // Read all genotype data
    let mut genotype_data = Vec::new();
    let mut sample_ids = Vec::new();
    
    for result in rdr.records() {
        let record = result.context("Failed to read genotype record")?;
        let sample_id = record[0].to_string();
        
        // Skip samples not in metadata
        if !sample_info.contains_key(&sample_id) {
            warn!("Sample {} not found in metadata, skipping", sample_id);
            continue;
        }
        
        let genotypes: Vec<f64> = record.iter().skip(1)
            .map(|x| x.parse::<f64>().unwrap_or(f64::NaN))
            .collect();
        
        genotype_data.push(genotypes);
        sample_ids.push(sample_id);
    }
    
    let total_samples = genotype_data.len();
    info!("Loaded genotypes for {} samples", total_samples);
    
    // Convert to ndarray
    let rows = genotype_data.len();
    let cols = genotype_data[0].len();
    let mut matrix = Array2::<f64>::zeros((rows, cols));
    
    for (i, row_data) in genotype_data.into_iter().enumerate() {
        for (j, val) in row_data.into_iter().enumerate() {
            matrix[[i, j]] = val;
        }
    }
    
    info!("Created genotype matrix: {} samples × {} variants", rows, cols);
    
    // Quality control filtering
    let (filtered_matrix, filtered_variants, qc_metrics) = 
        apply_quality_filters(&matrix, &variant_ids, &sample_ids, args)?;
    
    Ok((filtered_matrix, sample_ids, filtered_variants, qc_metrics))
}

fn apply_quality_filters(
    matrix: &Array2<f64>,
    variant_ids: &[String],
    sample_ids: &[String],
    args: &Args,
) -> Result<(Array2<f64>, Vec<String>, QualityMetrics)> {
    
    info!("Applying quality control filters...");
    
    let (n_samples, n_variants) = matrix.dim();
    let mut keep_variants = vec![true; n_variants];
    let mut keep_samples = vec![true; n_samples];
    
    // Calculate missing data rates
    let mut missing_per_variant = vec![0.0; n_variants];
    let mut missing_per_sample = vec![0.0; n_samples];
    
    // Calculate missing rates per variant
    for j in 0..n_variants {
        let missing_count = matrix.column(j).iter()
            .filter(|&&x| x.is_nan())
            .count();
        missing_per_variant[j] = missing_count as f64 / n_samples as f64;
        
        if missing_per_variant[j] > args.missing_threshold {
            keep_variants[j] = false;
        }
    }
    
    // Calculate missing rates per sample
    for i in 0..n_samples {
        let missing_count = matrix.row(i).iter()
            .filter(|&&x| x.is_nan())
            .count();
        missing_per_sample[i] = missing_count as f64 / n_variants as f64;
        
        if missing_per_sample[i] > args.sample_missing_threshold {
            keep_samples[i] = false;
        }
    }
    
    // Filter by minor allele frequency
    for j in 0..n_variants {
        if !keep_variants[j] { continue; }
        
        let col = matrix.column(j);
        let valid_genotypes: Vec<f64> = col.iter()
            .filter(|&&x| !x.is_nan())
            .copied()
            .collect();
        
        if valid_genotypes.is_empty() {
            keep_variants[j] = false;
            continue;
        }
        
        // Calculate allele frequency (assuming 0, 1, 2 encoding)
        let total_alleles = valid_genotypes.len() * 2;
        let alt_allele_count: f64 = valid_genotypes.iter().sum();
        let alt_freq = alt_allele_count / (total_alleles as f64);
        let maf = alt_freq.min(1.0 - alt_freq);
        
        if maf < args.maf_threshold {
            keep_variants[j] = false;
        }
    }
    
    // Create filtered matrix
    let variants_to_keep: Vec<usize> = keep_variants.iter()
        .enumerate()
        .filter(|(_, &keep)| keep)
        .map(|(i, _)| i)
        .collect();
    
    let samples_to_keep: Vec<usize> = keep_samples.iter()
        .enumerate()
        .filter(|(_, &keep)| keep)
        .map(|(i, _)| i)
        .collect();
    
    let filtered_variants: Vec<String> = variants_to_keep.iter()
        .map(|&i| variant_ids[i].clone())
        .collect();
    
    let filtered_matrix = matrix.select(Axis(0), &samples_to_keep)
        .select(Axis(1), &variants_to_keep);
    
    // Replace NaN with mean imputation
    let mut imputed_matrix = filtered_matrix.clone();
    for j in 0..imputed_matrix.ncols() {
        let mut col = imputed_matrix.column_mut(j);
        let valid_values: Vec<f64> = col.iter()
            .filter(|&&x| !x.is_nan())
            .copied()
            .collect();
        
        if !valid_values.is_empty() {
            let mean = valid_values.iter().sum::<f64>() / valid_values.len() as f64;
            for val in col.iter_mut() {
                if val.is_nan() {
                    *val = mean;
                }
            }
        }
    }
    
    let qc_metrics = QualityMetrics {
        total_variants: n_variants,
        variants_after_maf: keep_variants.iter().filter(|&&x| x).count(),
        variants_after_missing: variants_to_keep.len(),
        total_samples: n_samples,
        samples_after_qc: samples_to_keep.len(),
        mean_missing_per_sample: missing_per_sample.iter().sum::<f64>() / n_samples as f64,
        mean_missing_per_variant: missing_per_variant.iter().sum::<f64>() / n_variants as f64,
    };
    
    info!("QC Summary:");
    info!("  Variants: {} → {} ({:.1}% retained)", 
          qc_metrics.total_variants, qc_metrics.variants_after_missing,
          100.0 * qc_metrics.variants_after_missing as f64 / qc_metrics.total_variants as f64);
    info!("  Samples: {} → {} ({:.1}% retained)", 
          qc_metrics.total_samples, qc_metrics.samples_after_qc,
          100.0 * qc_metrics.samples_after_qc as f64 / qc_metrics.total_samples as f64);
    
    Ok((imputed_matrix, filtered_variants, qc_metrics))
}

#[derive(Debug)]
struct PCAResults {
    eigenvalues: Array1<f64>,
    eigenvectors: Array2<f64>,
    explained_variance_ratio: Array1<f64>,
    cumulative_variance: Array1<f64>,
}

fn perform_pca(matrix: &Array2<f64>, num_pcs: usize) -> Result<PCAResults> {
    let (n_samples, n_variants) = matrix.dim();
    
    // Center the matrix
    let mean_vec = matrix.mean_axis(Axis(0))
        .ok_or_else(|| AnalysisError::MatrixOperation("Failed to calculate mean".to_string()))?;
    
    let mut centered_matrix = matrix.clone();
    for mut row in centered_matrix.rows_mut() {
        row -= &mean_vec;
    }
    
    // Standardize (optional, but often recommended for genetic data)
    for j in 0..n_variants {
        let mut col = centered_matrix.column_mut(j);
        let std_dev = col.mapv(|x| x * x).sum().sqrt() / ((n_samples - 1) as f64).sqrt();
        if std_dev > 1e-8 {
            col /= std_dev;
        }
    }
    
    info!("Computing SVD for PCA...");
    
    // For large matrices, use covariance matrix approach
    let cov_matrix = if n_samples < n_variants {
        // Use sample covariance matrix (n_samples x n_samples)
        centered_matrix.dot(&centered_matrix.t()) / (n_variants as f64 - 1.0)
    } else {
        // Use variant covariance matrix approach
        let gram = centered_matrix.t().dot(&centered_matrix) / (n_samples as f64 - 1.0);
        gram
    };
    
    // Eigendecomposition
    let eigen_result = ndarray_linalg::eigh::eigh::Eigh::eigh(cov_matrix, ndarray_linalg::UPLO::Upper)
        .map_err(|e| AnalysisError::MatrixOperation(format!("Eigendecomposition failed: {:?}", e)))?;
    
    let (mut eigenvalues, mut eigenvectors) = (eigen_result.0, eigen_result.1);
    
    // Sort by eigenvalues (descending)
    let mut sorted_indices: Vec<usize> = (0..eigenvalues.len()).collect();
    sorted_indices.sort_by(|&i, &j| eigenvalues[j].partial_cmp(&eigenvalues[i]).unwrap());
    
    let n_components = num_pcs.min(eigenvalues.len());
    let selected_indices = &sorted_indices[..n_components];
    
    let selected_eigenvalues = Array1::from_vec(
        selected_indices.iter().map(|&i| eigenvalues[i]).collect()
    );
    
    let selected_eigenvectors = Array2::from_shape_vec(
        (eigenvectors.nrows(), n_components),
        selected_indices.iter()
            .flat_map(|&i| eigenvectors.column(i).iter().copied())
            .collect()
    ).map_err(|e| AnalysisError::MatrixOperation(format!("Failed to create eigenvector matrix: {}", e)))?;
    
    // Calculate explained variance
    let total_variance: f64 = eigenvalues.iter().filter(|&&x| x > 0.0).sum();
    let explained_variance_ratio = &selected_eigenvalues / total_variance;
    
    let mut cumulative_variance = explained_variance_ratio.clone();
    for i in 1..cumulative_variance.len() {
        cumulative_variance[i] += cumulative_variance[i - 1];
    }
    
    info!("PCA completed. Explained variance by first {} components:", n_components);
    for i in 0..n_components.min(5) {
        info!("  PC{}: {:.2}% ({:.2}% cumulative)", 
              i + 1, 
              explained_variance_ratio[i] * 100.0,
              cumulative_variance[i] * 100.0);
    }
    
    Ok(PCAResults {
        eigenvalues: selected_eigenvalues,
        eigenvectors: selected_eigenvectors,
        explained_variance_ratio,
        cumulative_variance,
    })
}

fn save_qc_report(output_dir: &str, metrics: &QualityMetrics) -> Result<()> {
    let report_path = format!("{}/qc_report.txt", output_dir);
    let mut file = File::create(&report_path)?;
    
    writeln!(file, "Population Genomics QC Report")?;
    writeln!(file, "==============================")?;
    writeln!(file, "")?;
    writeln!(file, "Sample Statistics:")?;
    writeln!(file, "  Total samples: {}", metrics.total_samples)?;
    writeln!(file, "  Samples after QC: {}", metrics.samples_after_qc)?;
    writeln!(file, "  Retention rate: {:.2}%", 
             100.0 * metrics.samples_after_qc as f64 / metrics.total_samples as f64)?;
    writeln!(file, "  Mean missing rate per sample: {:.4}", metrics.mean_missing_per_sample)?;
    writeln!(file, "")?;
    writeln!(file, "Variant Statistics:")?;
    writeln!(file, "  Total variants: {}", metrics.total_variants)?;
    writeln!(file, "  Variants after MAF filter: {}", metrics.variants_after_maf)?;
    writeln!(file, "  Variants after missing filter: {}", metrics.variants_after_missing)?;
    writeln!(file, "  Final retention rate: {:.2}%", 
             100.0 * metrics.variants_after_missing as f64 / metrics.total_variants as f64)?;
    writeln!(file, "  Mean missing rate per variant: {:.4}", metrics.mean_missing_per_variant)?;
    
    info!("QC report saved to {}", report_path);
    Ok(())
}

fn save_pca_results(
    output_dir: &str,
    pca_results: &PCAResults,
    sample_ids: &[String],
    sample_info: &HashMap<String, SampleInfo>,
) -> Result<()> {
    // Save eigenvalues and explained variance
    let eigenvalues_path = format!("{}/pca_eigenvalues.csv", output_dir);
    let mut wtr = Writer::from_path(&eigenvalues_path)?;
    wtr.write_record(&["PC", "Eigenvalue", "ExplainedVariance", "CumulativeVariance"])?;
    
    for (i, (&eigenval, (&exp_var, &cum_var))) in pca_results.eigenvalues.iter()
        .zip(pca_results.explained_variance_ratio.iter()
        .zip(pca_results.cumulative_variance.iter()))
        .enumerate() {
        wtr.write_record(&[
            &format!("PC{}", i + 1),
            &format!("{:.6}", eigenval),
            &format!("{:.6}", exp_var),
            &format!("{:.6}", cum_var),
        ])?;
    }
    wtr.flush()?;
    
    // Save PC scores for samples
    let scores_path = format!("{}/pca_scores.csv", output_dir);
    let mut wtr = Writer::from_path(&scores_path)?;
    
    let mut header = vec!["SAMPLE_ID".to_string(), "POPULATION".to_string()];
    for i in 0..pca_results.eigenvectors.ncols() {
        header.push(format!("PC{}", i + 1));
    }
    wtr.write_record(&header)?;
    
    for (sample_idx, sample_id) in sample_ids.iter().enumerate() {
        let mut record = vec![sample_id.clone()];
        
        // Add population info
        let population = sample_info.get(sample_id)
            .map(|info| info.population.clone())
            .unwrap_or_else(|| "Unknown".to_string());
        record.push(population);
        
        // Add PC scores
        for pc_idx in 0..pca_results.eigenvectors.ncols() {
            let score = pca_results.eigenvectors[[sample_idx, pc_idx]];
            record.push(format!("{:.6}", score));
        }
        
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    
    info!("PCA results saved to {} and {}", eigenvalues_path, scores_path);
    Ok(())
}

fn calculate_population_structure(
    matrix: &Array2<f64>,
    sample_ids: &[String],
    sample_info: &HashMap<String, SampleInfo>,
    output_dir: &str,
) -> Result<()> {
    // Group samples by population
    let mut pop_groups: HashMap<String, Vec<usize>> = HashMap::new();
    
    for (sample_idx, sample_id) in sample_ids.iter().enumerate() {
        if let Some(info) = sample_info.get(sample_id) {
            pop_groups.entry(info.population.clone())
                .or_insert_with(Vec::new)
                .push(sample_idx);
        }
    }
    
    let allele_freq_path = format!("{}/population_allele_frequencies.csv", output_dir);
    let mut wtr = Writer::from_path(&allele_freq_path)?;
    
    let mut header = vec!["VARIANT_ID".to_string()];
    let populations: Vec<String> = pop_groups.keys().cloned().collect();
    for pop in &populations {
        header.push(format!("{}_FREQ", pop));
        header.push(format!("{}_COUNT", pop));
    }
    wtr.write_record(&header)?;
    
    for variant_idx in 0..matrix.ncols() {
        let mut record = vec![format!("SNP_{:06}", variant_idx + 1)];
        
        for pop in &populations {
            if let Some(sample_indices) = pop_groups.get(pop) {
                let genotypes: Vec<f64> = sample_indices.iter()
                    .map(|&idx| matrix[[idx, variant_idx]])
                    .filter(|&x| !x.is_nan())
                    .collect();
                
                if !genotypes.is_empty() {
                    let total_alleles = genotypes.len() * 2;
                    let alt_allele_count: f64 = genotypes.iter().sum();
                    let alt_freq = alt_allele_count / (total_alleles as f64);
                    
                    record.push(format!("{:.6}", alt_freq));
                    record.push(genotypes.len().to_string());
                } else {
                    record.push("NA".to_string());
                    record.push("0".to_string());
                }
            }
        }
        
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    
    // Calculate Fst between populations
    calculate_fst(&matrix, sample_ids, sample_info, output_dir, &pop_groups)?;
    
    info!("Population structure analysis saved to {}", allele_freq_path);
    Ok(())
}

fn calculate_fst(
    matrix: &Array2<f64>,
    _sample_ids: &[String],
    _sample_info: &HashMap<String, SampleInfo>,
    output_dir: &str,
    pop_groups: &HashMap<String, Vec<usize>>,
) -> Result<()> {
    let fst_path = format!("{}/fst_matrix.csv", output_dir);
    let mut wtr = Writer::from_path(&fst_path)?;
    
    let populations: Vec<String> = pop_groups.keys().cloned().collect();
    let mut header = vec!["Population".to_string()];
    header.extend(populations.iter().cloned());
    wtr.write_record(&header)?;
    
    for pop1 in &populations {
        let mut record = vec![pop1.clone()];
        
        for pop2 in &populations {
            if pop1 == pop2 {
                record.push("0.0".to_string());
                continue;
            }
            
            let indices1 = &pop_groups[pop1];
            let indices2 = &pop_groups[pop2];
            
            let mut total_fst = 0.0;
            let mut valid_variants = 0;
            
            for variant_idx in 0..matrix.ncols() {
                if let Some(fst) = calculate_variant_fst(matrix, variant_idx, indices1, indices2) {
                    total_fst += fst;
                    valid_variants += 1;
                }
            }
            
            let mean_fst = if valid_variants > 0 {
                total_fst / valid_variants as f64
            } else {
                0.0
            };
            
            record.push(format!("{:.6}", mean_fst));
        }
        
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    
    info!("Fst matrix saved to {}", fst_path);
    Ok(())
}

fn calculate_variant_fst(
    matrix: &Array2<f64>,
    variant_idx: usize,
    indices1: &[usize],
    indices2: &[usize],
) -> Option<f64> {
    // Calculate allele frequencies for both populations
    let genotypes1: Vec<f64> = indices1.iter()
        .map(|&idx| matrix[[idx, variant_idx]])
        .filter(|&x| !x.is_nan())
        .collect();
    
    let genotypes2: Vec<f64> = indices2.iter()
        .map(|&idx| matrix[[idx, variant_idx]])
        .filter(|&x| !x.is_nan())
        .collect();
    
    if genotypes1.is_empty() || genotypes2.is_empty() {
        return None;
    }
    
    let n1 = genotypes1.len();
    let n2 = genotypes2.len();
    
    let p1 = genotypes1.iter().sum::<f64>() / (2.0 * n1 as f64);
    let p2 = genotypes2.iter().sum::<f64>() / (2.0 * n2 as f64);
    
    let p_total = (genotypes1.iter().sum::<f64>() + genotypes2.iter().sum::<f64>()) 
        / (2.0 * (n1 + n2) as f64);
    
    // Calculate heterozygosity
    let h1 = 2.0 * p1 * (1.0 - p1);
    let h2 = 2.0 * p2 * (1.0 - p2);
    let ht = 2.0 * p_total * (1.0 - p_total);
    
    let hs = (n1 as f64 * h1 + n2 as f64 * h2) / (n1 + n2) as f64;
    
    if ht > 0.0 {
        Some((ht - hs) / ht)
    } else {
        None
    }
}

fn perform_hwe_tests(
    matrix: &Array2<f64>,
    variant_ids: &[String],
    output_dir: &str,
) -> Result<()> {
    let hwe_path = format!("{}/hwe_tests.csv", output_dir);
    let mut wtr = Writer::from_path(&hwe_path)?;
    wtr.write_record(&["VARIANT_ID", "OBSERVED_HET", "EXPECTED_HET", "HWE_PVALUE", "SIGNIFICANT"])?;
    
    for (variant_idx, variant_id) in variant_ids.iter().enumerate() {
        let genotypes: Vec<i32> = matrix.column(variant_idx)
            .iter()
            .filter_map(|&x| if x.is_nan() { None } else { Some(x as i32) })
            .collect();
        
        if genotypes.is_empty() {
            continue;
        }
        
        let (obs_het, exp_het, p_value) = hardy_weinberg_test(&genotypes);
        let significant = if p_value < 0.001 { "***" } 
                         else if p_value < 0.01 { "**" } 
                         else if p_value < 0.05 { "*" } 
                         else { "" };
        
        wtr.write_record(&[
            variant_id,
            &format!("{:.4}", obs_het),
            &format!("{:.4}", exp_het),
            &format!("{:.2e}", p_value),
            significant,
        ])?;
    }
    wtr.flush()?;
    
    info!("Hardy-Weinberg equilibrium tests saved to {}", hwe_path);
    Ok(())
}

fn hardy_weinberg_test(genotypes: &[i32]) -> (f64, f64, f64) {
    let mut counts = [0; 3]; // AA, Aa, aa
    
    for &gt in genotypes {
        if gt >= 0 && gt <= 2 {
            counts[gt as usize] += 1;
        }
    }
    
    let total = counts.iter().sum::<i32>() as f64;
    if total == 0.0 {
        return (0.0, 0.0, 1.0);
    }
    
    let n_aa = counts[0] as f64;
    let n_ab = counts[1] as f64;
    let n_bb = counts[2] as f64;
    
    // Calculate allele frequency
    let p = (2.0 * n_aa + n_ab) / (2.0 * total);
    let q = 1.0 - p;
    
    // Expected genotype frequencies under HWE
    let exp_aa = total * p * p;
    let exp_ab = total * 2.0 * p * q;
    let exp_bb = total * q * q;
    
    // Observed and expected heterozygosity
    let obs_het = n_ab / total;
    let exp_het = 2.0 * p * q;
    
    // Chi-square test
    let chi_square = if exp_aa > 0.0 && exp_ab > 0.0 && exp_bb > 0.0 {
        (n_aa - exp_aa).powi(2) / exp_aa +
        (n_ab - exp_ab).powi(2) / exp_ab +
        (n_bb - exp_bb).powi(2) / exp_bb
    } else {
        0.0
    };
    
    // P-value (using chi-square distribution with 1 df)
    let p_value = if chi_square > 0.0 {
        use statrs::distribution::{ChiSquared, ContinuousCDF};
        let chi_sq_dist = ChiSquared::new(1.0).unwrap();
        1.0 - chi_sq_dist.cdf(chi_square)
    } else {
        1.0
    };
    
    (obs_het, exp_het, p_value)
}