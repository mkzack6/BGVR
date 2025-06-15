#!/usr/bin/env python3
"""
Generate realistic genomic data for testing the enhanced genomic algorithms pipeline.

This script creates population genetics datasets with realistic patterns including:
- Hardy-Weinberg equilibrium deviations
- Linkage disequilibrium blocks
- Population structure
- Realistic allele frequency distributions
- Missing data patterns

Output format: CSV with individuals as rows and variants as columns
Values: 0 (homozygous reference), 1 (heterozygous), 2 (homozygous alternate), 9 (missing)
"""

import random
import numpy as np
import argparse
from pathlib import Path
import json
from datetime import datetime

class PopulationGenomicsSimulator:
    def __init__(self, seed=42):
        """Initialize the population genetics simulator."""
        random.seed(seed)
        np.random.seed(seed)
        self.populations = []
        self.ld_blocks = []
        
    def generate_allele_frequencies(self, n_variants, distribution='realistic'):
        """Generate realistic allele frequency distributions."""
        if distribution == 'uniform':
            return np.random.uniform(0.05, 0.95, n_variants)
        elif distribution == 'realistic':
            # Realistic MAF distribution: many rare variants, fewer common
            frequencies = []
            for _ in range(n_variants):
                if random.random() < 0.4:  # 40% rare variants
                    freq = np.random.beta(0.5, 5)  # Skewed toward low frequencies
                elif random.random() < 0.8:  # 40% intermediate frequency
                    freq = np.random.beta(2, 2)    # More uniform
                else:  # 20% common variants
                    freq = np.random.beta(2, 1)    # Skewed toward high frequencies
                frequencies.append(max(0.01, min(0.99, freq)))
            return np.array(frequencies)
        else:
            raise ValueError(f"Unknown distribution: {distribution}")
    
    def create_ld_blocks(self, n_variants, block_size_range=(10, 100), n_blocks=None):
        """Create linkage disequilibrium blocks."""
        if n_blocks is None:
            n_blocks = max(1, n_variants // 50)  # Roughly one block per 50 variants
        
        blocks = []
        remaining_variants = list(range(n_variants))
        
        for _ in range(n_blocks):
            if len(remaining_variants) < 5:
                break
                
            block_size = random.randint(*block_size_range)
            block_size = min(block_size, len(remaining_variants))
            
            # Select consecutive variants for the block
            start_idx = random.randint(0, len(remaining_variants) - block_size)
            block_variants = remaining_variants[start_idx:start_idx + block_size]
            
            # Remove selected variants from remaining
            remaining_variants = [v for v in remaining_variants if v not in block_variants]
            
            # Create LD structure within block
            r_squared_matrix = np.random.exponential(0.3, (block_size, block_size))
            r_squared_matrix = np.minimum(r_squared_matrix, 1.0)
            
            # Make symmetric and set diagonal to 1
            for i in range(block_size):
                for j in range(i, block_size):
                    if i == j:
                        r_squared_matrix[i, j] = 1.0
                    else:
                        val = r_squared_matrix[i, j]
                        r_squared_matrix[i, j] = val
                        r_squared_matrix[j, i] = val
            
            blocks.append({
                'variants': block_variants,
                'r_squared_matrix': r_squared_matrix
            })
        
        return blocks
    
    def generate_haplotypes_with_ld(self, n_individuals, allele_frequencies, ld_blocks):
        """Generate haplotypes with realistic LD patterns."""
        n_variants = len(allele_frequencies)
        haplotypes = np.zeros((n_individuals * 2, n_variants), dtype=np.uint8)
        
        # Track which variants are in LD blocks
        block_variants = set()
        for block in ld_blocks:
            block_variants.update(block['variants'])
        
        # Generate independent variants first
        independent_variants = [i for i in range(n_variants) if i not in block_variants]
        for var_idx in independent_variants:
            freq = allele_frequencies[var_idx]
            haplotypes[:, var_idx] = np.random.binomial(1, freq, n_individuals * 2)
        
        # Generate LD block variants
        for block in ld_blocks:
            variants = block['variants']
            if len(variants) == 0:
                continue
                
            block_size = len(variants)
            
            # Generate correlated haplotypes for this block
            # Start with the first variant
            first_var = variants[0]
            freq = allele_frequencies[first_var]
            haplotypes[:, first_var] = np.random.binomial(1, freq, n_individuals * 2)
            
            # Generate subsequent variants with LD
            for i, var_idx in enumerate(variants[1:], 1):
                freq = allele_frequencies[var_idx]
                
                # Calculate correlation with previous variants in block
                correlated_prob = np.zeros(n_individuals * 2)
                for hap_idx in range(n_individuals * 2):
                    # Base probability
                    prob = freq
                    
                    # Adjust based on correlation with other variants in block
                    for j in range(i):
                        prev_var = variants[j]
                        if j < len(block['r_squared_matrix']) and i < len(block['r_squared_matrix'][j]):
                            r_sq = block['r_squared_matrix'][j][i]
                            # Simple LD model: if previous variant is 1, increase probability
                            if haplotypes[hap_idx, prev_var] == 1:
                                prob = prob + r_sq * (1 - prob) * 0.5
                            else:
                                prob = prob * (1 - r_sq * 0.5)
                    
                    correlated_prob[hap_idx] = max(0.01, min(0.99, prob))
                
                # Generate alleles with correlated probabilities
                for hap_idx in range(n_individuals * 2):
                    haplotypes[hap_idx, var_idx] = np.random.binomial(1, correlated_prob[hap_idx])
        
        return haplotypes
    
    def combine_haplotypes_to_genotypes(self, haplotypes):
        """Combine paired haplotypes into diploid genotypes."""
        n_total_haps, n_variants = haplotypes.shape
        n_individuals = n_total_haps // 2
        
        genotypes = np.zeros((n_individuals, n_variants), dtype=np.uint8)
        
        for ind_idx in range(n_individuals):
            hap1_idx = ind_idx * 2
            hap2_idx = ind_idx * 2 + 1
            genotypes[ind_idx, :] = haplotypes[hap1_idx, :] + haplotypes[hap2_idx, :]
        
        return genotypes
    
    def add_missing_data(self, genotypes, missing_rate=0.02):
        """Add realistic missing data patterns."""
        n_individuals, n_variants = genotypes.shape
        
        # Individual-level missingness (some individuals have more missing data)
        individual_missing_rates = np.random.beta(1, 20, n_individuals) * missing_rate * 10
        individual_missing_rates = np.clip(individual_missing_rates, 0, 0.2)
        
        # Variant-level missingness (some variants harder to call)
        variant_missing_rates = np.random.beta(1, 10, n_variants) * missing_rate * 5
        variant_missing_rates = np.clip(variant_missing_rates, 0, 0.1)
        
        # Apply missing data
        for i in range(n_individuals):
            for j in range(n_variants):
                total_missing_prob = individual_missing_rates[i] + variant_missing_rates[j]
                if random.random() < total_missing_prob:
                    genotypes[i, j] = 9  # Missing data marker
        
        return genotypes
    
    def add_population_structure(self, genotypes, n_populations=3):
        """Add population structure by modifying allele frequencies."""
        n_individuals, n_variants = genotypes.shape
        
        # Assign individuals to populations
        pop_assignments = np.random.multinomial(1, [1/n_populations] * n_populations, n_individuals)
        pop_labels = np.argmax(pop_assignments, axis=1)
        
        # Create population-specific allele frequency shifts
        population_effects = np.random.normal(0, 0.1, (n_populations, n_variants))
        
        # Apply population effects
        for i in range(n_individuals):
            pop = pop_labels[i]
            for j in range(n_variants):
                if genotypes[i, j] != 9:  # Don't modify missing data
                    # Apply small population-specific effects
                    effect = population_effects[pop, j]
                    if random.random() < abs(effect):
                        if effect > 0 and genotypes[i, j] < 2:
                            genotypes[i, j] += 1
                        elif effect < 0 and genotypes[i, j] > 0:
                            genotypes[i, j] -= 1
        
        return genotypes, pop_labels
    
    def generate_dataset(self, 
                        n_individuals=1000, 
                        n_variants=5000,
                        missing_rate=0.02,
                        n_populations=1,
                        ld_blocks=True,
                        maf_distribution='realistic'):
        """Generate a complete genomic dataset."""
        
        print(f"Generating dataset: {n_individuals} individuals × {n_variants} variants")
        
        # Generate allele frequencies
        allele_frequencies = self.generate_allele_frequencies(n_variants, maf_distribution)
        print(f"Generated allele frequencies (mean MAF: {np.mean(np.minimum(allele_frequencies, 1-allele_frequencies)):.4f})")
        
        # Create LD blocks if requested
        ld_block_structures = []
        if ld_blocks:
            ld_block_structures = self.create_ld_blocks(n_variants)
            print(f"Created {len(ld_block_structures)} LD blocks")
        
        # Generate haplotypes
        haplotypes = self.generate_haplotypes_with_ld(n_individuals, allele_frequencies, ld_block_structures)
        print("Generated haplotypes with LD structure")
        
        # Combine to diploid genotypes
        genotypes = self.combine_haplotypes_to_genotypes(haplotypes)
        print("Combined haplotypes to diploid genotypes")
        
        # Add missing data
        genotypes = self.add_missing_data(genotypes, missing_rate)
        missing_count = np.sum(genotypes == 9)
        total_genotypes = n_individuals * n_variants
        actual_missing_rate = missing_count / total_genotypes
        print(f"Added missing data ({actual_missing_rate:.4f} missing rate)")
        
        # Add population structure if requested
        pop_labels = None
        if n_populations > 1:
            genotypes, pop_labels = self.add_population_structure(genotypes, n_populations)
            print(f"Added population structure ({n_populations} populations)")
        
        return {
            'genotypes': genotypes,
            'allele_frequencies': allele_frequencies,
            'ld_blocks': ld_block_structures,
            'population_labels': pop_labels,
            'metadata': {
                'n_individuals': n_individuals,
                'n_variants': n_variants,
                'missing_rate': actual_missing_rate,
                'n_populations': n_populations,
                'mean_maf': float(np.mean(np.minimum(allele_frequencies, 1-allele_frequencies))),
                'n_ld_blocks': len(ld_block_structures)
            }
        }

def write_csv_genotypes(genotypes, output_file):
    """Write genotypes to CSV format."""
    print(f"Writing genotypes to {output_file}")
    
    n_individuals, n_variants = genotypes.shape
    
    with open(output_file, 'w') as f:
        # Write header comment
        f.write("# Genomic genotype data in CSV format\n")
        f.write(f"# Individuals: {n_individuals}, Variants: {n_variants}\n")
        f.write("# Format: 0=homozygous ref, 1=heterozygous, 2=homozygous alt, 9=missing\n")
        f.write("# Rows=individuals, Columns=variants\n")
        
        # Write genotype data
        for i in range(n_individuals):
            row = ','.join(str(genotypes[i, j]) for j in range(n_variants))
            f.write(row + '\n')

def write_metadata(dataset, output_file):
    """Write dataset metadata to JSON file."""
    metadata = dataset['metadata'].copy()
    metadata['generation_date'] = datetime.now().isoformat()
    metadata['ld_blocks_summary'] = [
        {
            'block_id': i,
            'n_variants': len(block['variants']),
            'variant_range': f"{min(block['variants'])}-{max(block['variants'])}"
        }
        for i, block in enumerate(dataset['ld_blocks'][:10])  # First 10 blocks only
    ]
    
    with open(output_file, 'w') as f:
        json.dump(metadata, f, indent=2)

def main():
    parser = argparse.ArgumentParser(
        description="Generate realistic genomic datasets for testing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-o", "--output",
        default="genotypes.csv",
        help="Output CSV file name"
    )
    
    parser.add_argument(
        "-i", "--individuals",
        type=int,
        default=1000,
        help="Number of individuals"
    )
    
    parser.add_argument(
        "-v", "--variants",
        type=int,
        default=5000,
        help="Number of variants"
    )
    
    parser.add_argument(
        "-m", "--missing-rate",
        type=float,
        default=0.02,
        help="Missing data rate"
    )
    
    parser.add_argument(
        "-p", "--populations",
        type=int,
        default=1,
        help="Number of populations"
    )
    
    parser.add_argument(
        "--no-ld",
        action="store_true",
        help="Disable LD block structure"
    )
    
    parser.add_argument(
        "--maf-distribution",
        choices=['uniform', 'realistic'],
        default='realistic',
        help="Minor allele frequency distribution"
    )
    
    parser.add_argument(
        "-s", "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility"
    )
    
    parser.add_argument(
        "--metadata",
        help="Output metadata JSON file"
    )
    
    # Preset configurations
    preset_group = parser.add_mutually_exclusive_group()
    preset_group.add_argument(
        "--small",
        action="store_true",
        help="Small dataset: 100 individuals, 1000 variants"
    )
    
    preset_group.add_argument(
        "--medium",
        action="store_true",
        help="Medium dataset: 1000 individuals, 5000 variants"
    )
    
    preset_group.add_argument(
        "--large",
        action="store_true",
        help="Large dataset: 5000 individuals, 20000 variants"
    )
    
    args = parser.parse_args()
    
    # Apply presets
    if args.small:
        args.individuals = 100
        args.variants = 1000
    elif args.medium:
        args.individuals = 1000
        args.variants = 5000
    elif args.large:
        args.individuals = 5000
        args.variants = 20000
    
    # Generate dataset
    simulator = PopulationGenomicsSimulator(seed=args.seed)
    
    dataset = simulator.generate_dataset(
        n_individuals=args.individuals,
        n_variants=args.variants,
        missing_rate=args.missing_rate,
        n_populations=args.populations,
        ld_blocks=not args.no_ld,
        maf_distribution=args.maf_distribution
    )
    
    # Write outputs
    write_csv_genotypes(dataset['genotypes'], args.output)
    
    if args.metadata:
        write_metadata(dataset, args.metadata)
    
    # Print summary
    print("\n=== Dataset Generation Summary ===")
    print(f"Output file: {args.output}")
    print(f"Individuals: {dataset['metadata']['n_individuals']}")
    print(f"Variants: {dataset['metadata']['n_variants']}")
    print(f"Missing rate: {dataset['metadata']['missing_rate']:.4f}")
    print(f"Mean MAF: {dataset['metadata']['mean_maf']:.4f}")
    print(f"LD blocks: {dataset['metadata']['n_ld_blocks']}")
    print(f"Populations: {dataset['metadata']['n_populations']}")
    print(f"File size: {Path(args.output).stat().st_size / 1024 / 1024:.2f} MB")

if __name__ == "__main__":
    main()