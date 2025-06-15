#!/usr/bin/env python3
"""
Generate sample CSV files for population genomics pipeline
This script creates genotypes.csv and sample_metadata.csv files
that can be used immediately with the pipeline.
"""

import pandas as pd
import numpy as np
import argparse
import os
from pathlib import Path

def generate_population_data(n_individuals=1000, n_snps=10000, n_populations=3, output_dir="data"):
    """
    Generate synthetic population genomics data
    
    Parameters:
    - n_individuals: Number of individuals to simulate
    - n_snps: Number of SNPs to generate
    - n_populations: Number of populations
    - output_dir: Output directory for files
    """
    
    print(f"🧬 Generating population genomics data...")
    print(f"   Individuals: {n_individuals}")
    print(f"   SNPs: {n_snps}")
    print(f"   Populations: {n_populations}")
    print(f"   Output: {output_dir}/")
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Generate population-specific parameters
    populations = []
    for i in range(n_populations):
        pop_name = f"POP{i+1}"
        populations.append({
            'name': pop_name,
            'drift': 0.1 + i * 0.05,  # Increasing drift
            'size': n_individuals // n_populations
        })
    
    # Generate base allele frequencies
    base_allele_freqs = np.random.beta(2, 2, n_snps)  # Beta distribution for realistic frequencies
    base_allele_freqs = np.clip(base_allele_freqs, 0.05, 0.95)  # Keep in reasonable range
    
    # Generate population-specific allele frequencies
    pop_allele_freqs = {}
    for pop in populations:
        # Add population-specific drift
        drift_effect = np.random.normal(0, pop['drift'], n_snps)
        pop_freqs = base_allele_freqs + drift_effect
        pop_freqs = np.clip(pop_freqs, 0.01, 0.99)
        pop_allele_freqs[pop['name']] = pop_freqs
    
    # Generate sample metadata
    print("Generating sample metadata...")
    sample_data = []
    current_pop_idx = 0
    individuals_in_current_pop = 0
    
    for i in range(n_individuals):
        # Switch to next population if current one is full
        if (individuals_in_current_pop >= populations[current_pop_idx]['size'] and 
            current_pop_idx < n_populations - 1):
            current_pop_idx += 1
            individuals_in_current_pop = 0
        
        sample_id = f"IND_{i+1:06d}"
        population = populations[current_pop_idx]['name']
        sex = np.random.choice(['M', 'F'])
        age = np.random.randint(18, 80)
        # Generate correlated phenotype (slightly different between populations)
        phenotype_base = 50 + current_pop_idx * 10
        phenotype = np.random.normal(phenotype_base, 15)
        
        sample_data.append({
            'SAMPLE_ID': sample_id,
            'POPULATION': population,
            'SEX': sex,
            'AGE': age,
            'PHENOTYPE': round(phenotype, 2)
        })
        
        individuals_in_current_pop += 1
    
    # Create sample metadata DataFrame
    sample_df = pd.DataFrame(sample_data)
    
    # Generate genotype data
    print("Generating genotype matrix...")
    genotype_data = []
    
    for i, sample in enumerate(sample_data):
        if i % 100 == 0:
            print(f"   Generated genotypes for {i} individuals...")
        
        population = sample['POPULATION']
        sample_id = sample['SAMPLE_ID']
        
        # Generate genotypes for this individual
        genotypes = [sample_id]  # Start with sample ID
        
        for snp_idx in range(n_snps):
            allele_freq = pop_allele_freqs[population][snp_idx]
            
            # Generate diploid genotype (0, 1, or 2 copies of reference allele)
            # Using binomial distribution
            genotype = np.random.binomial(2, allele_freq)
            genotypes.append(genotype)
        
        genotype_data.append(genotypes)
    
    # Create genotype DataFrame
    print("Creating genotype matrix...")
    snp_columns = [f"SNP_{i+1:06d}" for i in range(n_snps)]
    genotype_columns = ['SAMPLE_ID'] + snp_columns
    genotype_df = pd.DataFrame(genotype_data, columns=genotype_columns)
    
    # Generate variant metadata
    print("Generating variant metadata...")
    variant_data = []
    chromosomes = list(range(1, 23))  # Chromosomes 1-22
    alleles = ['A', 'T', 'G', 'C']
    
    for i in range(n_snps):
        snp_id = f"SNP_{i+1:06d}"
        chromosome = np.random.choice(chromosomes)
        position = np.random.randint(1000000, 250000000)
        ref_allele = np.random.choice(alleles)
        alt_allele = np.random.choice([a for a in alleles if a != ref_allele])
        gene = f"GENE_{np.random.randint(1, 5000):04d}"
        
        variant_data.append({
            'SNP_ID': snp_id,
            'CHROMOSOME': chromosome,
            'POSITION': position,
            'REF_ALLELE': ref_allele,
            'ALT_ALLELE': alt_allele,
            'GENE': gene
        })
    
    variant_df = pd.DataFrame(variant_data)
    
    # Save files
    print("Saving files...")
    
    # Save genotypes
    genotype_file = os.path.join(output_dir, "genotypes.csv")
    genotype_df.to_csv(genotype_file, index=False)
    print(f"✅ Genotypes saved to: {genotype_file}")
    
    # Save sample metadata
    metadata_file = os.path.join(output_dir, "sample_metadata.csv")
    sample_df.to_csv(metadata_file, index=False)
    print(f"✅ Sample metadata saved to: {metadata_file}")
    
    # Save variant metadata
    variant_file = os.path.join(output_dir, "variant_metadata.csv")
    variant_df.to_csv(variant_file, index=False)
    print(f"✅ Variant metadata saved to: {variant_file}")
    
    # Generate summary statistics
    print("\n📊 Data Summary:")
    print(f"   Total samples: {len(sample_df)}")
    print(f"   Total SNPs: {n_snps}")
    print("   Population distribution:")
    pop_counts = sample_df['POPULATION'].value_counts()
    for pop, count in pop_counts.items():
        print(f"     {pop}: {count} individuals")
    
    print(f"   Age range: {sample_df['AGE'].min()}-{sample_df['AGE'].max()}")
    print(f"   Sex distribution: {dict(sample_df['SEX'].value_counts())}")
    print(f"   Phenotype range: {sample_df['PHENOTYPE'].min():.2f}-{sample_df['PHENOTYPE'].max():.2f}")
    
    # Calculate some basic statistics
    print("\n🔍 Genotype Statistics:")
    
    # Sample a few SNPs for statistics (to avoid memory issues with large datasets)
    sample_snps = min(1000, n_snps)
    sampled_genotypes = genotype_df.iloc[:, 1:sample_snps+1].values
    
    # Calculate missing rate (assuming -1 or NaN would be missing)
    missing_rate = np.sum(pd.isna(sampled_genotypes)) / sampled_genotypes.size
    print(f"   Missing rate: {missing_rate:.4f}")
    
    # Calculate mean genotype value
    mean_genotype = np.nanmean(sampled_genotypes)
    print(f"   Mean genotype value: {mean_genotype:.3f}")
    
    # Calculate allele frequency distribution
    allele_freqs = np.nanmean(sampled_genotypes, axis=0) / 2.0
    print(f"   Mean allele frequency: {np.mean(allele_freqs):.3f}")
    print(f"   Allele frequency range: {np.min(allele_freqs):.3f}-{np.max(allele_freqs):.3f}")
    
    print(f"\n🎉 Data generation completed successfully!")
    print(f"Files are ready to use with the population genomics pipeline.")
    
    return genotype_file, metadata_file, variant_file

def main():
    parser = argparse.ArgumentParser(description="Generate synthetic population genomics data")
    parser.add_argument("--individuals", "-i", type=int, default=1000,
                       help="Number of individuals (default: 1000)")
    parser.add_argument("--snps", "-s", type=int, default=10000,
                       help="Number of SNPs (default: 10000)")
    parser.add_argument("--populations", "-p", type=int, default=3,
                       help="Number of populations (default: 3)")
    parser.add_argument("--output", "-o", type=str, default="data",
                       help="Output directory (default: data)")
    parser.add_argument("--small", action="store_true",
                       help="Generate small test dataset (100 individuals, 1000 SNPs)")
    parser.add_argument("--large", action="store_true",
                       help="Generate large dataset (5000 individuals, 50000 SNPs)")
    
    args = parser.parse_args()
    
    # Predefined dataset sizes
    if args.small:
        individuals, snps, populations = 100, 1000, 2
        print("🧪 Generating SMALL test dataset...")
    elif args.large:
        individuals, snps, populations = 5000, 50000, 5
        print("🚀 Generating LARGE dataset...")
    else:
        individuals, snps, populations = args.individuals, args.snps, args.populations
    
    try:
        generate_population_data(
            n_individuals=individuals,
            n_snps=snps,
            n_populations=populations,
            output_dir=args.output
        )
    except Exception as e:
        print(f"❌ Error generating data: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())