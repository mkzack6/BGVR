# 6.4. Parallel and Distributed Processing of HTS Data
## 6.4.1: Parallel BAM Read Counter with Rust and Nextflow

This project provides a Rust program (`bam_counter`) integrated with a Nextflow pipeline to perform parallelized read counting across multiple genomic regions in BAM files. Designed for educational purposes, it demonstrates high-performance bioinformatics workflows using Rust for efficient read counting and Nextflow for scalable pipeline orchestration. The program processes BAM files, counts total and mapped reads in specified genomic regions, and is suitable for analyzing large-scale genomic datasets in high-performance computing (HPC) environments.

## Breakdown of the Project

### Dependencies

The project consists of a Rust program (`bam_counter`) and a Nextflow pipeline (`main.nf`). Dependencies are defined in `Cargo.toml` for the Rust program and system-level tools for the pipeline:

#### Rust Dependencies (Cargo.toml):
```toml
[dependencies]
clap = { version = "4.0", features = ["derive"] }
anyhow = "1.0"
log = "0.4"
rayon = "1.7"
rust-htslib = "0.49.0"
```

- **anyhow**: Robust error handling for BAM file I/O and record processing
- **clap**: Command-line parsing for specifying BAM files and regions
- **log**: Structured logging for debugging and metrics
- **rayon**: Parallel processing for efficient read counting across regions
- **rust-htslib**: Core BAM parsing and manipulation

#### System Dependencies:
- **Nextflow**: Orchestrates the pipeline, managing parallel execution across genomic regions
- **samtools**: Used for viewing and validating BAM files
- **Rust/Cargo**: Required to build and run the Rust program

### Implementation Modules

#### Rust Program (src/main.rs):

**Data Structures:**
- `ReadCounts`: Tracks total and mapped read counts
- `Cli`: Command-line interface parameters parsed via clap

**Read Counting Logic:**
- Uses rust-htslib for BAM file reading
- Counts total reads and mapped reads in specified regions
- Processes regions in parallel using rayon

**Input/Output:**
- Input: BAM file and list of regions (e.g., 'chr1:1-1000000')
- Output: Read count statistics for each region

**Main Function:**
- Parses command-line arguments
- Processes regions in parallel
- Aggregates and reports results

#### Nextflow Pipeline (main.nf):

**Workflow:**
- Defines input channels for sample processing
- Implements three main processes: alignmentOrFetch, variantCalling, and mergeVcfs
- Handles both mock and real data processing modes

**Processes:**
- `alignmentOrFetch`: Generates or fetches BAM files
- `variantCalling`: Processes BAM files using the Rust tool
- `mergeVcfs`: Aggregates results

**Input/Output:**
- Input: Sample list (samples.txt)
- Output: Processed VCF files in results directory

### Program Execution

#### Prerequisites

1. Install Rust/Cargo:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

2. Install Nextflow:
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

3. Install samtools:
```bash
sudo apt update
sudo apt install samtools
```

#### Setup

1. Clone or create project directory:
```bash
mkdir bam_counter
cd bam_counter
```

2. Build Rust Program:
```bash
cargo build --release
```

3. Verify the binary:
```bash
ls target/release/bam_counter
```

#### Run Program

1. Create sample list:
```bash
echo -e "sample1\nsample2\nsample3" > samples.txt
```

2. Run Nextflow Pipeline:
```bash
./nextflow run main.nf
```

### Example Output

```
N E X T F L O W  ~  version 24.10.5
Launching `main.nf` [elegant_noyce] DSL2 - revision: 3b457e4cd7

executor >  local (7)
[6c/5350f4] alignmentOrFetch (3) [100%] 3 of 3 ✔
[3d/ee2120] variantCalling (2)   [100%] 3 of 3 ✔
[22/d7f871] mergeVcfs            [100%] 1 of 1 ✔
```

### Potential Errors and Solutions

1. **BAM File Not Found:**
   - Error: `Failed to open indexed BAM file`
   - Solution: Ensure BAM file exists and is indexed (`samtools index file.bam`)

2. **Invalid Region Format:**
   - Error: `Failed to fetch region`
   - Solution: Use correct format (e.g., 'chr1:1-1000000')

3. **Missing BAM Index:**
   - Error: `No index file found`
   - Solution: Create index with `samtools index`

### Why is this Project Important?

1. **Bioinformatics Workflow:** Demonstrates parallel processing of BAM files, essential for genomic analysis pipelines
2. **High-Performance Computing:** Leverages Rust's speed and Nextflow's scalability
3. **Educational Value:** Teaches modern bioinformatics tools and parallel processing
4. **Scalability:** Suitable for HPC environments and large datasets
5. **Error Handling:** Demonstrates robust error handling in bioinformatics pipelines

### Next Steps / Improvements

1. **Enhanced Metrics:**
   - Add more detailed read statistics (e.g., mapping quality distribution)
   - Generate QC reports

2. **Pipeline Extensions:**
   - Add support for paired-end read analysis
   - Implement coverage calculation
   - Add variant calling capabilities

3. **Performance Optimization:**
   - Implement streaming processing for large files
   - Add memory usage optimization options

4. **Visualization:**
   - Add coverage plot generation
   - Create interactive HTML reports

5. **Additional Features:**
   - Support for CRAM format
   - Integration with cloud storage
   - Multi-sample comparison tools 
