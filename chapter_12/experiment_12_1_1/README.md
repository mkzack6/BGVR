# Population Genomics Analysis Pipeline

A comprehensive Rust-based pipeline for population genomics analysis, orchestrated with Nextflow. This pipeline provides end-to-end functionality for synthetic data generation, quality control, principal component analysis, population structure analysis, and Hardy-Weinberg equilibrium testing.

![Pipeline Overview](https://img.shields.io/badge/Pipeline-Population%20Genomics-blue) ![Language](https://img.shields.io/badge/Language-Rust-orange) ![Workflow](https://img.shields.io/badge/Workflow-Nextflow-green) ![Platform](https://img.shields.io/badge/Platform-WSL%2FLinux-lightgrey)

## 🧬 Features

- *Synthetic Data Generation*: Create realistic population genomics datasets with configurable population structure
- *Quality Control*: Comprehensive SNP and sample filtering based on missing data rates and minor allele frequency
- *Principal Component Analysis*: Standardized PCA with eigenvalue decomposition for population structure visualization
- *Population Structure Analysis*: Population-specific allele frequency analysis and pairwise Fst calculations
- *Hardy-Weinberg Testing*: Statistical testing for equilibrium deviations across variants
- *Reproducible Workflows*: Nextflow orchestration ensuring scalability and reproducibility
- *Comprehensive Reporting*: Automated HTML and text report generation

## 🚀 Quick Start

### Prerequisites
- Ubuntu/Debian Linux (WSL2 recommended for Windows users)
- Internet connection for dependency installation
- Minimum 4GB RAM (8GB+ recommended for large datasets)
- 2GB free disk space

### 1. Automated Setup
bash
# Download and run the setup script
curl -sSL https://raw.githubusercontent.com/your-repo/setup.sh | bash


*OR manually:*
bash
# Clone/download the setup script
wget https://raw.githubusercontent.com/your-repo/setup.sh
chmod +x setup.sh
./setup.sh


### 2. Navigate to Project
bash
cd population_genomics_pipeline


### 3. Copy Source Code
*Important*: Copy the source code from the artifacts into the appropriate files:

bash
# Copy the Rust source files (paste artifact content):
nano src/bin/generate_data.rs       # Copy generate_data.rs artifact
nano src/bin/population_analysis.rs # Copy population_analysis.rs artifact
nano main.nf                        # Copy main.nf artifact


### 4. Build and Test
bash
# Compile Rust binaries
cargo build --release

# Run quick test with small dataset
./scripts/test_pipeline.sh


### 5. Run Full Analysis
bash
# Standard pipeline run
./scripts/run_pipeline.sh

# Or with custom parameters
./scripts/run_pipeline.sh --individuals 2000 --snps 20000 --populations 5


## 📁 Directory Structure


population_genomics_pipeline/
├── Cargo.toml                    # Rust dependencies and project configuration
├── main.nf                       # Nextflow pipeline definition
├── nextflow.config               # Pipeline execution configuration
├── README.md                     # This documentation file
├── LICENSE                       # MIT license
├── src/
│   └── bin/
│       ├── generate_data.rs      # Synthetic data generator
│       └── population_analysis.rs # Main analysis pipeline
├── scripts/                      # Helper scripts
│   ├── run_pipeline.sh           # Main pipeline execution script
│   ├── test_pipeline.sh          # Quick test runner
│   └── generate_data_only.sh     # Data generation only
├── data/                         # Input data directory
│   ├── genotypes.csv             # Generated genotype matrix
│   ├── sample_metadata.csv       # Sample information
│   └── variant_metadata.csv      # Variant annotations
├── results/                      # Analysis output directory
│   ├── qc/                       # Quality control results
│   ├── pca/                      # PCA analysis results
│   ├── population_structure/     # Population genetics results
│   ├── hwe/                      # Hardy-Weinberg test results
│   ├── final_report.html         # Comprehensive HTML report
│   └── summary_statistics.txt    # Text summary
├── target/                       # Rust compilation artifacts
└── work/                         # Nextflow execution cache


## 🔧 Pipeline Components

### 1. Data Generation (generate_data.rs)
Creates synthetic population genomics datasets with realistic characteristics:

*Features:*
- Configurable number of individuals, SNPs, and populations
- Population-specific allele frequency drift
- Diploid genotype encoding (0/1/2 format)
- Comprehensive metadata generation
- Random phenotype simulation

*Usage:*
bash
./target/release/generate_data \
    --individuals 1000 \
    --snps 10000 \
    --populations 3 \
    --output data/


### 2. Population Analysis (population_analysis.rs)
Comprehensive genomic analysis pipeline:

*Quality Control:*
- Minor allele frequency (MAF) filtering
- Missing data rate filtering (per variant and per sample)
- Sample quality assessment
- Comprehensive QC reporting

*Principal Component Analysis:*
- Matrix standardization and centering
- Singular value decomposition (SVD)
- Explained variance calculation
- Population structure visualization data

*Population Structure Analysis:*
- Population-specific allele frequency calculation
- Pairwise Fst estimation using Weir & Cockerham method
- Population differentiation metrics

*Hardy-Weinberg Equilibrium Testing:*
- Chi-square goodness-of-fit tests
- Multiple testing considerations
- Deviation significance assessment

*Usage:*
bash
./target/release/population_analysis \
    --genotypes data/genotypes.csv \
    --metadata data/sample_metadata.csv \
    --output results/ \
    --maf-threshold 0.05 \
    --missing-threshold 0.1 \
    --num-pcs 10 \
    --hwe-test \
    --population-structure


### 3. Nextflow Orchestration (main.nf)
Scalable workflow management:

*Process Flow:*
1. *GENERATE_DATA*: Synthetic dataset creation
2. *QUALITY_CONTROL*: Data filtering and QC assessment
3. *PCA_ANALYSIS*: Principal component computation
4. *POPULATION_STRUCTURE*: Population genetics analysis (optional)
5. *HWE_TESTING*: Hardy-Weinberg equilibrium testing (optional)
6. *GENERATE_REPORT*: Comprehensive report compilation

*Resource Management:*
- Configurable CPU and memory allocation
- Process-specific resource requirements
- Execution monitoring and reporting

## ⚙ Configuration Parameters

### Data Generation Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| --individuals | 1000 | Number of individuals to simulate |
| --snps | 10000 | Number of SNPs to generate |
| --populations | 3 | Number of distinct populations |
| --output | "data" | Output directory for generated data |

### Analysis Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| --genotypes | "data/genotypes.csv" | Input genotype file path |
| --metadata | "data/sample_metadata.csv" | Sample metadata file path |
| --output | "results" | Analysis output directory |
| --maf_threshold | 0.05 | Minor allele frequency threshold |
| --missing_threshold | 0.1 | Missing data threshold per variant |
| --sample_missing_threshold | 0.1 | Missing data threshold per sample |
| --num_pcs | 10 | Number of principal components to compute |
| --hwe_test | true | Enable Hardy-Weinberg equilibrium testing |
| --population_structure | true | Enable population structure analysis |

### Nextflow Parameters
All analysis parameters can be passed to Nextflow with double dashes:
bash
nextflow run main.nf --individuals 2000 --maf_threshold 0.01 --num_pcs 20


## 📊 Output Files

### Quality Control (results/qc/)
- **qc_report.txt**: Comprehensive quality control summary
- **qc_summary.log**: Detailed filtering process log

### PCA Analysis (results/pca/)
- **pca_eigenvalues.csv**: Eigenvalues, explained variance, and cumulative variance
- **pca_scores.csv**: Sample coordinates in principal component space
- **pca_analysis.log**: PCA computation log

### Population Structure (results/population_structure/)
- **population_allele_frequencies.csv**: Population-specific allele frequencies
- **fst_matrix.csv**: Pairwise Fst values between populations
- **population_structure.log**: Analysis process log

### Hardy-Weinberg Testing (results/hwe/)
- **hwe_tests.csv**: HWE test results with p-values and significance flags
- **hwe_analysis.log**: Testing process log

### Summary Reports (results/)
- **final_report.html**: Interactive HTML report with all results
- **summary_statistics.txt**: Consolidated text summary
- **pipeline_report.html**: Nextflow execution report
- **timeline.html**: Pipeline execution timeline
- **dag.svg**: Directed acyclic graph of pipeline processes

## 🚀 Usage Examples

### Basic Pipeline Execution
bash
# Standard run with default parameters
nextflow run main.nf

# Quick test with small dataset
nextflow run main.nf --individuals 100 --snps 1000 --populations 2


### Large-Scale Analysis
bash
# Large cohort analysis
nextflow run main.nf \
    --individuals 5000 \
    --snps 50000 \
    --populations 8 \
    --maf_threshold 0.01 \
    --num_pcs 20


### Focused Analysis
bash
# Only PCA analysis, skip population structure
nextflow run main.nf \
    --population_structure false \
    --hwe_test false \
    --num_pcs 15


### Using Helper Scripts
bash
# Quick test
./scripts/test_pipeline.sh

# Custom run with scripts
./scripts/run_pipeline.sh --individuals 3000 --snps 30000

# Generate data only
./scripts/generate_data_only.sh 2000 20000 4


### Individual Component Usage
bash
# Generate data only
cargo run --bin generate_data -- \
    --individuals 1000 --snps 10000 --populations 3

# Analysis only (with existing data)
cargo run --bin population_analysis -- \
    --genotypes data/genotypes.csv \
    --metadata data/sample_metadata.csv \
    --hwe-test --population-structure


## 🔧 Advanced Configuration

### Resource Management
Edit nextflow.config for custom resource allocation:

groovy
process {
    executor = 'local'
    memory = '8 GB'
    cpus = 4
    
    withName:GENERATE_DATA {
        memory = '4 GB'
        cpus = 2
    }
    
    withName:PCA_ANALYSIS {
        memory = '16 GB'
        cpus = 8
        time = '4h'
    }
}


### Cluster Execution
For cluster environments, add cluster profiles:

groovy
profiles {
    cluster {
        process.executor = 'slurm'
        process.queue = 'normal'
        process.memory = '16 GB'
        process.time = '2h'
    }
}


Run with: nextflow run main.nf -profile cluster

### Container Support
Enable Docker or Singularity in nextflow.config:

groovy
docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
}

process.container = 'your-registry/population-genomics:latest'


## 🐛 Troubleshooting

### Common Issues

#### Compilation Errors
bash
# Install system dependencies
sudo apt update
sudo apt install build-essential pkg-config libssl-dev libblas-dev liblapack-dev

# Update Rust toolchain
rustup update stable


#### Memory Issues
bash
# Check available memory
free -h

# Reduce dataset size for testing
./scripts/test_pipeline.sh

# Increase swap space if needed
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile


#### Nextflow Issues
bash
# Check Java version (requires Java 11+)
java -version

# Clean Nextflow cache
rm -rf work/ .nextflow*

# Update Nextflow
nextflow self-update


### Debug Mode
Enable verbose logging:

bash
# Rust applications
RUST_LOG=debug ./target/release/population_analysis [options]

# Nextflow pipeline
nextflow run main.nf -with-trace -with-report -with-timeline


### Performance Monitoring
bash
# Monitor system resources
htop

# Check disk space
df -h

# Monitor pipeline progress
tail -f .nextflow.log


## 📈 Performance Characteristics

### Computational Complexity
- *Data Generation*: O(individuals × SNPs)
- *Quality Control*: O(individuals × SNPs)
- *PCA*: O(min(individuals, SNPs)³)
- *Population Structure*: O(populations² × SNPs)
- *HWE Testing*: O(SNPs)

### Memory Requirements
| Dataset Size | Recommended RAM |
|--------------|-----------------|
| 1K ind × 10K SNPs | 4 GB |
| 5K ind × 50K SNPs | 16 GB |
| 10K ind × 100K SNPs | 32 GB |
| 50K ind × 500K SNPs | 128 GB |

### Runtime Estimates
| Dataset Size | Approximate Runtime |
|--------------|-------------------|
| 1K ind × 10K SNPs | 2-5 minutes |
| 5K ind × 50K SNPs | 15-30 minutes |
| 10K ind × 100K SNPs | 1-2 hours |
| 50K ind × 500K SNPs | 4-8 hours |

Times vary significantly based on hardware specifications

## 🤝 Contributing

### Development Setup
bash
# Clone the repository
git clone https://github.com/your-repo/population-genomics-pipeline
cd population-genomics-pipeline

# Install development dependencies
cargo install cargo-watch cargo-audit

# Run tests
cargo test

# Format code
cargo fmt

# Lint code
cargo clippy


### Code Style
- Follow Rust standard formatting (cargo fmt)
- Use meaningful variable names
- Add comprehensive documentation
- Include unit tests for new functions
- Follow error handling best practices

### Submitting Changes
1. Fork the repository
2. Create a feature branch
3. Make your changes with tests
4. Ensure all tests pass
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- *Rust Community*: For excellent libraries (ndarray, tokio, clap)
- *Nextflow Team*: For robust workflow orchestration
- *Population Genetics Community*: For algorithmic foundations
- *Contributors*: All contributors to this project

## 📚 References

- Di Tommaso, P. et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- Weir, B. S., & Cockerham, C. C. (1984). Estimating F-statistics for the analysis of population structure. Evolution, 38(6), 1358-1370.
- Patterson, N., Price, A. L., & Reich, D. (2006). Population structure and eigenanalysis. PLoS Genetics, 2(12), e190.

## 📞 Support

For questions, issues, or feature requests:

1. *Check Documentation*: Review this README and inline code documentation
2. *Search Issues*: Look for existing GitHub issues
3. *Create Issue*: Submit a detailed issue report
4. *Community Forum*: Join discussions in our community forum

---

*Happy analyzing! 🧬🔬*
