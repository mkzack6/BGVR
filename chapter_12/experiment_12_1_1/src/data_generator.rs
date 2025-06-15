use csv::Writer;
use rand::prelude::*;
use rand_distr::{Binomial, Normal};
use std::fs::File;
use std::io::Write;
use clap::Parser;
use anyhow::Result;

/// Generate synthetic population genomics data
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Number of individuals
    #[arg(short, long, default_value_t = 1000)]
    individuals: usize,
    
    /// Number of SNPs
    #[arg(short, long, default_value_t = 10000)]
    snps: usize,
    
    /// Number of populations
    #[arg(short, long, default_value_t = 3)]
    populations: usize,
    
    /// Output directory
    #[arg(short, long, default_value = "data")]
    output: String,
}

#[derive(Debug)]
struct Population {
    name: String,
    size: usize,
    allele_frequencies: Vec<f64>,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();
    
    // Create output directory
    std::fs::create_dir_all(&args.output)?;
    
    log::info!("Generating synthetic population genomics data...");
    log::info!("Individuals: {}, SNPs: {}, Populations: {}", 
               args.individuals, args.snps, args.populations);
    
    // Generate populations with different allele frequencies
    let populations = generate_populations(args.populations, args.snps)?;
    
    // Generate genotype data
    generate_genotype_data(&args, &populations)?;
    
    // Generate sample metadata
    generate_sample_metadata(&args, &populations)?;
    
    // Generate variant metadata
    generate_variant_metadata(&args)?;
    
    log::info!("Data generation complete!");
    Ok(())
}

fn generate_populations(num_pops: usize, num_snps: usize) -> Result<Vec<Population>> {
    let mut rng = thread_rng();
    let mut populations = Vec::new();
    
    for i in 0..num_pops {
        let pop_name = format!("POP{}", i + 1);
        
        // Generate allele frequencies with population-specific drift
        let mut allele_freqs = Vec::with_capacity(num_snps);
        
        for _ in 0..num_snps {
            // Base frequency with some population-specific drift
            let base_freq = rng.gen_range(0.05..0.95);
            let drift = Normal::new(0.0, 0.1 + i as f64 * 0.05).unwrap();
            let pop_freq = (base_freq + drift.sample(&mut rng)).clamp(0.01, 0.99);
            allele_freqs.push(pop_freq);
        }
        
        populations.push(Population {
            name: pop_name,
            size: 0, // Will be filled later
            allele_frequencies: allele_freqs,
        });
    }
    
    Ok(populations)
}

fn generate_genotype_data(args: &Args, populations: &[Population]) -> Result<()> {
    let mut rng = thread_rng();
    let genotype_file = format!("{}/genotypes.csv", args.output);
    let mut wtr = Writer::from_path(&genotype_file)?;
    
    // Create header
    let mut header = vec!["SAMPLE_ID".to_string()];
    for i in 0..args.snps {
        header.push(format!("SNP_{:06}", i + 1));
    }
    wtr.write_record(&header)?;
    
    // Assign individuals to populations
    let individuals_per_pop = args.individuals / args.populations;
    let mut current_pop = 0;
    
    for ind_idx in 0..args.individuals {
        // Switch population every individuals_per_pop individuals
        if ind_idx > 0 && ind_idx % individuals_per_pop == 0 && current_pop < args.populations - 1 {
            current_pop += 1;
        }
        
        let sample_id = format!("IND_{:06}", ind_idx + 1);
        let mut genotypes = vec![sample_id];
        
        // Generate genotypes for this individual
        for snp_idx in 0..args.snps {
            let allele_freq = populations[current_pop].allele_frequencies[snp_idx];
            
            // Generate diploid genotype (0, 1, or 2 copies of reference allele)
            let binomial = Binomial::new(2, allele_freq).unwrap();
            let genotype = binomial.sample(&mut rng);
            genotypes.push(genotype.to_string());
        }
        
        wtr.write_record(&genotypes)?;
        
        if ind_idx % 100 == 0 {
            log::info!("Generated genotypes for {} individuals", ind_idx + 1);
        }
    }
    
    wtr.flush()?;
    log::info!("Genotype data written to {}", genotype_file);
    Ok(())
}

fn generate_sample_metadata(args: &Args, populations: &[Population]) -> Result<()> {
    let metadata_file = format!("{}/sample_metadata.csv", args.output);
    let mut wtr = Writer::from_path(&metadata_file)?;
    
    // Header
    wtr.write_record(&["SAMPLE_ID", "POPULATION", "SEX", "AGE", "PHENOTYPE"])?;
    
    let mut rng = thread_rng();
    let individuals_per_pop = args.individuals / args.populations;
    let mut current_pop = 0;
    
    for ind_idx in 0..args.individuals {
        if ind_idx > 0 && ind_idx % individuals_per_pop == 0 && current_pop < args.populations - 1 {
            current_pop += 1;
        }
        
        let sample_id = format!("IND_{:06}", ind_idx + 1);
        let population = &populations[current_pop].name;
        let sex = if rng.gen_bool(0.5) { "M" } else { "F" };
        let age = rng.gen_range(18..80);
        let phenotype = rng.gen_range(0.0..100.0);
        
        wtr.write_record(&[
            &sample_id,
            population,
            sex,
            &age.to_string(),
            &format!("{:.2}", phenotype),
        ])?;
    }
    
    wtr.flush()?;
    log::info!("Sample metadata written to {}", metadata_file);
    Ok(())
}

fn generate_variant_metadata(args: &Args) -> Result<()> {
    let variant_file = format!("{}/variant_metadata.csv", args.output);
    let mut wtr = Writer::from_path(&variant_file)?;
    
    // Header
    wtr.write_record(&["SNP_ID", "CHROMOSOME", "POSITION", "REF_ALLELE", "ALT_ALLELE", "GENE"])?;
    
    let mut rng = thread_rng();
    let chromosomes: Vec<&str> = (1..=22).map(|i| Box::leak(i.to_string().into_boxed_str())).collect();
    let alleles = ["A", "T", "G", "C"];
    
    for snp_idx in 0..args.snps {
        let snp_id = format!("SNP_{:06}", snp_idx + 1);
        let chromosome = chromosomes.choose(&mut rng).unwrap();
        let position = rng.gen_range(1000000..250000000);
        let ref_allele = alleles.choose(&mut rng).unwrap();
        let mut alt_allele = alleles.choose(&mut rng).unwrap();
        while alt_allele == ref_allele {
            alt_allele = alleles.choose(&mut rng).unwrap();
        }
        let gene = format!("GENE_{:04}", rng.gen_range(1..5000));
        
        wtr.write_record(&[
            &snp_id,
            chromosome,
            &position.to_string(),
            ref_allele,
            alt_allele,
            &gene,
        ])?;
    }
    
    wtr.flush()?;
    log::info!("Variant metadata written to {}", variant_file);
    Ok(())
}