#!/bin/bash

# Quick CSV Generator for Population Genomics Pipeline
# ===================================================
# This script generates the required CSV files without needing Rust compilation

set -e

# Default parameters
INDIVIDUALS=1000
SNPS=10000
POPULATIONS=3
OUTPUT_DIR="data"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --individuals|-i)
            INDIVIDUALS="$2"
            shift 2
            ;;
        --snps|-s)
            SNPS="$2"
            shift 2
            ;;
        --populations|-p)
            POPULATIONS="$2"
            shift 2
            ;;
        --output|-o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --small)
            INDIVIDUALS=100
            SNPS=1000
            POPULATIONS=2
            print_info "Using SMALL dataset parameters"
            shift
            ;;
        --test)
            INDIVIDUALS=50
            SNPS=500
            POPULATIONS=2
            print_info "Using TEST dataset parameters"
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --individuals, -i N    Number of individuals (default: $INDIVIDUALS)"
            echo "  --snps, -s N          Number of SNPs (default: $SNPS)"
            echo "  --populations, -p N   Number of populations (default: $POPULATIONS)"
            echo "  --output, -o DIR      Output directory (default: $OUTPUT_DIR)"
            echo "  --small               Quick small dataset (100 ind, 1000 SNPs)"
            echo "  --test                Tiny test dataset (50 ind, 500 SNPs)"
            echo "  --help, -h            Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo "🧬 Population Genomics CSV Generator"
echo "===================================="
echo "Parameters:"
echo "  Individuals: $INDIVIDUALS"
echo "  SNPs: $SNPS"
echo "  Populations: $POPULATIONS"
echo "  Output directory: $OUTPUT_DIR"
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    print_warning "Python3 not found. Trying python..."
    if ! command -v python &> /dev/null; then
        echo "❌ Python not found. Please install Python to use this script."
        echo "Alternative: Use the Rust data generator after compiling the project."
        exit 1
    else
        PYTHON_CMD="python"
    fi
else
    PYTHON_CMD="python3"
fi

# Check if required Python packages are available
print_step "Checking Python dependencies..."
$PYTHON_CMD -c "import pandas, numpy" 2>/dev/null || {
    print_warning "Required Python packages not found. Installing..."
    if command -v pip3 &> /dev/null; then
        pip3 install pandas numpy
    elif command -v pip &> /dev/null; then
        pip install pandas numpy
    else
        echo "❌ pip not found. Please install pandas and numpy manually:"
        echo "   sudo apt install python3-pandas python3-numpy"
        echo "   # OR"
        echo "   pip install pandas numpy"
        exit 1
    fi
}

# Create output directory
mkdir -p "$OUTPUT_DIR"

print_step "Generating CSV files..."

# Create the Python script inline and run it
$PYTHON_CMD << EOF
import pandas as pd
import numpy as np
import os

# Set parameters
n_individuals = $INDIVIDUALS
n_snps = $SNPS
n_populations = $POPULATIONS
output_dir = "$OUTPUT_DIR"

print(f"Generating {n_individuals} individuals, {n_snps} SNPs, {n_populations} populations...")

# Set random seed for reproducibility
np.random.seed(42)

# Generate population info
pop_names = [f"POP{i+1}" for i in range(n_populations)]
individuals_per_pop = n_individuals // n_populations

print("Generating sample metadata...")

# Generate sample metadata
sample_data = []
for i in range(n_individuals):
    pop_idx = i // individuals_per_pop
    if pop_idx >= n_populations:  # Handle remainder
        pop_idx = n_populations - 1
    
    sample_data.append({
        'SAMPLE_ID': f"IND_{i+1:06d}",
        'POPULATION': pop_names[pop_idx],
        'SEX': np.random.choice(['M', 'F']),
        'AGE': np.random.randint(18, 80),
        'PHENOTYPE': round(np.random.normal(50 + pop_idx * 10, 15), 2)
    })

sample_df = pd.DataFrame(sample_data)

print("Generating genotype data...")

# Generate population-specific allele frequencies
pop_allele_freqs = {}
base_freqs = np.random.beta(2, 2, n_snps)  # Realistic allele frequencies
base_freqs = np.clip(base_freqs, 0.05, 0.95)

for i, pop in enumerate(pop_names):
    # Add population-specific drift
    drift = np.random.normal(0, 0.1 + i * 0.05, n_snps)
    pop_freqs = np.clip(base_freqs + drift, 0.01, 0.99)
    pop_allele_freqs[pop] = pop_freqs

# Generate genotype matrix
genotype_data = []
for i, sample in enumerate(sample_data):
    if i % 100 == 0 and i > 0:
        print(f"  Generated {i} individuals...")
    
    population = sample['POPULATION']
    sample_id = sample['SAMPLE_ID']
    
    # Generate genotypes using population-specific frequencies
    genotypes = [sample_id]
    for snp_idx in range(n_snps):
        allele_freq = pop_allele_freqs[population][snp_idx]
        genotype = np.random.binomial(2, allele_freq)
        genotypes.append(genotype)
    
    genotype_data.append(genotypes)

# Create genotype DataFrame
snp_columns = [f"SNP_{i+1:06d}" for i in range(n_snps)]
genotype_columns = ['SAMPLE_ID'] + snp_columns
genotype_df = pd.DataFrame(genotype_data, columns=genotype_columns)

print("Generating variant metadata...")

# Generate variant metadata
variant_data = []
chromosomes = list(range(1, 23))
alleles = ['A', 'T', 'G', 'C']

for i in range(n_snps):
    ref_allele = np.random.choice(alleles)
    alt_allele = np.random.choice([a for a in alleles if a != ref_allele])
    
    variant_data.append({
        'SNP_ID': f"SNP_{i+1:06d}",
        'CHROMOSOME': np.random.choice(chromosomes),
        'POSITION': np.random.randint(1000000, 250000000),
        'REF_ALLELE': ref_allele,
        'ALT_ALLELE': alt_allele,
        'GENE': f"GENE_{np.random.randint(1, 5000):04d}"
    })

variant_df = pd.DataFrame(variant_data)

print("Saving files...")

# Save files
genotype_df.to_csv(os.path.join(output_dir, "genotypes.csv"), index=False)
sample_df.to_csv(os.path.join(output_dir, "sample_metadata.csv"), index=False)
variant_df.to_csv(os.path.join(output_dir, "variant_metadata.csv"), index=False)

print("\\n📊 Data Summary:")
print(f"  Samples: {len(sample_df)}")
print(f"  SNPs: {n_snps}")
print(f"  Populations: {dict(sample_df['POPULATION'].value_counts())}")
print(f"  File size estimates:")
print(f"    genotypes.csv: ~{(n_individuals * n_snps * 2) // 1024 // 1024} MB")
print(f"    sample_metadata.csv: ~{len(sample_df) * 50 // 1024} KB")
print(f"    variant_metadata.csv: ~{n_snps * 50 // 1024} KB")

EOF

# Verify files were created
print_step "Verifying generated files..."

if [[ -f "$OUTPUT_DIR/genotypes.csv" ]]; then
    GENOTYPE_LINES=$(wc -l < "$OUTPUT_DIR/genotypes.csv")
    print_info "✅ genotypes.csv created ($((GENOTYPE_LINES-1)) samples)"
else
    echo "❌ Failed to create genotypes.csv"
    exit 1
fi

if [[ -f "$OUTPUT_DIR/sample_metadata.csv" ]]; then
    METADATA_LINES=$(wc -l < "$OUTPUT_DIR/sample_metadata.csv")
    print_info "✅ sample_metadata.csv created ($((METADATA_LINES-1)) samples)"
else
    echo "❌ Failed to create sample_metadata.csv"
    exit 1
fi

if [[ -f "$OUTPUT_DIR/variant_metadata.csv" ]]; then
    VARIANT_LINES=$(wc -l < "$OUTPUT_DIR/variant_metadata.csv")
    print_info "✅ variant_metadata.csv created ($((VARIANT_LINES-1)) variants)"
else
    echo "❌ Failed to create variant_metadata.csv"
    exit 1
fi

echo ""
print_info "🎉 CSV files generated successfully!"
echo ""
echo "Generated files:"
echo "  📁 $OUTPUT_DIR/"
echo "    📄 genotypes.csv ($(du -h "$OUTPUT_DIR/genotypes.csv" | cut -f1))"
echo "    📄 sample_metadata.csv ($(du -h "$OUTPUT_DIR/sample_metadata.csv" | cut -f1))"
echo "    📄 variant_metadata.csv ($(du -h "$OUTPUT_DIR/variant_metadata.csv" | cut -f1))"
echo ""
echo "Sample preview:"
echo "🔍 First 5 lines of sample_metadata.csv:"
head -5 "$OUTPUT_DIR/sample_metadata.csv"
echo ""
echo "🔍 First 3 columns of genotypes.csv (first 3 lines):"
cut -d',' -f1-3 "$OUTPUT_DIR/genotypes.csv" | head -3
echo ""
echo "Ready to use with the population genomics pipeline!"
echo "Next steps:"
echo "  1. cargo build --release  # (if using Rust pipeline)"
echo "  2. ./scripts/test_pipeline.sh  # (test with generated data)"
echo "  3. nextflow run main.nf  # (run full analysis)"