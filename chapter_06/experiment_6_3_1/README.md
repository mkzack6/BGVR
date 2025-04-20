# 6.3. Variant Call Format (VCF/BCF) Handling
## Experiment 6.3.1: Parallel BCF Filtering with Rust and Nextflow

This project provides a Rust program (`bcf_filter_tool`) integrated with a Nextflow pipeline to filter genomic variant data in BCF format. Designed for educational purposes, it demonstrates high-performance bioinformatics workflows using Rust for efficient variant filtering and Nextflow for scalable pipeline orchestration. The program processes BCF files, applies quality and depth filters, and outputs filtered BCF files for specified genomic chunks. It is suitable for analyzing large-scale genomic datasets in high-performance computing (HPC) environments, showcasing parallel processing and robust error handling.

### Breakdown of the Project

#### Dependencies
The project consists of a Rust program (`bcf_filter_tool`) and a Nextflow pipeline (`main.nf`). Dependencies are defined in `Cargo.toml` for the Rust program and system-level tools for the pipeline:

**Rust Dependencies (`Cargo.toml`)**:
```toml
[dependencies]
anyhow = "1.0"
clap = { version = "4.4", features = ["derive"] }
env_logger = "0.11.7"
log = "0.4"
rayon = "1.8"
rust-htslib = "0.49.0"
```
- `bio`: Provides BCF file parsing and manipulation for genomic variant data.
- `clap`: Supports command-line argument parsing for filter parameters (e.g., `--min-qual`, `--min-depth`).
- `anyhow`: Facilitates flexible error handling for robust execution.

**System Dependencies**:
- `Nextflow`: Orchestrates the pipeline, managing parallel execution across genomic chunks.
- `bcftools`: Used for viewing and validating BCF files.
- `Rust/Cargo`: Required to build and run the Rust program.

#### Feature Flags
The program has no feature flags, using a single workflow with configurable parameters:
- **Default**: Expects a BCF file (`wgs_cohort.bcf`) and processes it in chunks (e.g., `chr1:1-2000`, `chr2:1-3000`).
- **Data**: Produces filtered BCF files (e.g., `filtered_chr1_1_2000.bcf`) with variants meeting quality and depth thresholds.
- **Parameters**:
  - `--min-qual 30`: Filters variants with QUAL < 30.
  - `--min-depth 10`: Filters variants with average depth < 10 across samples.
  - `--chunk-size 1000`: Defines the genomic region size for processing.

The pipeline has one main process:
- `filterVariants`: Runs `bcf_filter_tool` on each BCF file chunk, producing a filtered BCF file.

#### Implementation Modules
**Rust Program (`src/main.rs`)**:
- **Data Structures**:
  - `BCF Record`: Represents a genomic variant with fields like `CHROM`, `POS`, `QUAL`, and `DP` (depth) from the `bio` crate.
  - `FilterParams`: Stores filtering parameters (`min_qual`, `min_depth`) parsed via `clap`.
- **Filtering Logic**:
  - Reads BCF files using `bio::io::bcf`.
  - Applies filters: retains variants with `QUAL >= 30` and average `DP >= 10` across samples.
  - Writes filtered variants to a new BCF file.
- **Input/Output**:
  - **Input**: BCF file (`wgs_cohort.bcf`) and chunk (e.g., `chr1:1-2000`).
  - **Output**: Filtered BCF file (e.g., `filtered_chr1_1_2000.bcf`).
- **Main Function**:
  - Parses command-line arguments, reads the input BCF, applies filters, and writes the output BCF.
  - Logs metrics (e.g., number of records processed/retained) via `RUST_LOG=info`.

**Nextflow Pipeline (`main.nf`)**:
- **Workflow**:
  - Defines input channels for the BCF file and chunks.
  - Combines channels to create tuples `[bcf_file, chunk]`.
  - Runs `filterVariants` process in parallel for each tuple.
- **Process**:
  - Executes `bcf_filter_tool` with specified parameters.
  - Publishes filtered BCF files to `/home/zack/experiment_6_3_1/bcf_filter/`.
- **Input/Output**:
  - **Input**: BCF file (`wgs_cohort.bcf`), chunk list (e.g., `chr1:1-2000`, `chr2:1-3000`).
  - **Output**: Filtered BCF files named `filtered_<chunk>.bcf` (e.g., `filtered_chr1_1_2000.bcf`).

### Program Execution
The pipeline processes a BCF file to filter variants in parallel across genomic chunks, producing filtered BCF files. Steps to execute:

#### Prerequisites
1. **Install Rust/Cargo**:
   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   source $HOME/.cargo/env
   ```
2. **Install Nextflow**:
   ```bash
   curl -s https://get.nextflow.io | bash
   sudo mv nextflow /usr/local/bin/
   ```
3. **Install bcftools**:
   ```bash
   sudo apt update
   sudo apt install bcftools
   ```
4. **Setup Project Directory**:
   Ensure the project directory (`~/experiment_6_3_1/bcf_filter/`) contains:
   - `Cargo.toml` and `src/main.rs` for `bcf_filter_tool`.
   - `main.nf` for the Nextflow pipeline.
   - `wgs_cohort.bcf` and `wgs_cohort.bcf.csi` (BCF file and index).

#### Setup
1. **Verify Input BCF**:
   Ensure `wgs_cohort.bcf` exists and is indexed:
   ```bash
   ls ~/experiment_6_3_1/bcf_filter/wgs_cohort.bcf ~/experiment_6_3_1/bcf_filter/wgs_cohort.bcf.csi
   bcftools view ~/experiment_6_3_1/bcf_filter/wgs_cohort.bcf
   ```
   If the index is missing:
   ```bash
   bcftools index ~/experiment_6_3_1/bcf_filter/wgs_cohort.bcf
   ```
2. **Build Rust Program**:
   ```bash
   cd ~/experiment_6_3_1/bcf_filter
   cargo build --release
   ```
   Verify binary:
   ```bash
   ls target/release/bcf_filter_tool
   ```

#### Run Program
1. **Navigate to Project Directory**:
   ```bash
   cd ~/experiment_6_3_1/bcf_filter
   ```
2. **Run Nextflow Pipeline**:
   ```bash
   export RUST_LOG=info
   nextflow run main.nf
   ```
   **Note**: Ensure `main.nf` is configured with the correct paths for `wgs_cohort.bcf` and `bcf_filter_tool`.

#### Example Output
```
N E X T F L O W  ~  version 24.10.5
Launching `main.nf` [some_name] DSL2 - revision: ...
executor >  local (2)
[xx/xxxxxx] filterVariants (1) | 1 of 2
[xx/xxxxxx] filterVariants (2) | 2 of 2 ✔
/home/zack/experiment_6_3_1/bcf_filter/filtered_chr1_1_2000.bcf
/home/zack/experiment_6_3_1/bcf_filter/filtered_chr2_1_3000.bcf
```
**Rust Logs**:
- For `chr1:1-2000`:
  ```
  [INFO] Starting parallel BCF filtering with parameters: FilterParams { min_qual: 30.0, min_depth: 10 }
  [INFO] Filtering completed. Processed 2 records in total. Retained 1 records after filtering.
  [INFO] Filtering completed successfully.
  ```
- For `chr2:1-3000`:
  ```
  [INFO] Starting parallel BCF filtering with parameters: FilterParams { min_qual: 30.0, min_depth: 10 }
  [INFO] Filtering completed. Processed 1 records in total. Retained 1 records after filtering.
  [INFO] Filtering completed successfully.
  ```

#### Generated Files
- `filtered_chr1_1_2000.bcf`: Filtered BCF file for `chr1:1-2000` (contains 1 record: `chr1:1000`).
- `filtered_chr2_1_3000.bcf`: Filtered BCF file for `chr2:1-3000` (contains 1 record: `chr2:3000`).
- Location: `/home/zack/experiment_6_3_1/bcf_filter/`.

### Potential Errors and Solutions
1. **Channel Mismatch in Nextflow**:
   - **Error**: `Process `filterVariants` declares 2 input channels but 1 were specified`.
   - **Cause**: Incorrect tuple handling in `main.nf`.
   - **Solution**: Ensure `main.nf` uses `tuple path(bcf_file), val(chunk)` and `map { bcf_file, chunk -> tuple(bcf_file, chunk) }`.
2. **BCF File Not Found**:
   - **Error**: `Failed to open file "wgs_cohort.bcf" : No such file or directory`.
   - **Cause**: Missing or incorrect path to `wgs_cohort.bcf`.
   - **Solution**: Verify file exists and update `params.input_bcf` in `main.nf`.
3. **Incorrect Filenames**:
   - **Error**: Outputs like `filtered_chr1_1-2000.bcf` instead of `filtered_chr1_1_2000.bcf`.
   - **Cause**: Missing hyphen replacement in `chunk` string.
   - **Solution**: Use `chunk.replace(':','_').replace('-','_')` in `main.nf`.
4. **Files Not Published**:
   - **Error**: Files remain in `work/` directory.
   - **Cause**: `publishDir` directive misconfigured or permissions issue.
   - **Solution**: Verify `publishDir` path and directory permissions:
     ```bash
     ls -ld ~/experiment_6_3_1/bcf_filter
     chmod 775 ~/experiment_6_3_1/bcf_filter
     ```
5. **Rust Binary Failure**:
   - **Error**: `Error processing chunk ...`.
   - **Cause**: Rust program fails to process BCF (e.g., invalid chunk or BCF format).
   - **Solution**: Test manually:
     ```bash
     ./target/release/bcf_filter_tool --input wgs_cohort.bcf --output test.bcf --min-qual 30 --min-depth 10 --chunk-size 1000
     ```

### Why is this Project Important?
- **Bioinformatics Workflow**: Demonstrates parallel filtering of genomic variants, essential for variant calling pipelines in genomics research.
- **High-Performance Computing**: Leverages Rust’s speed and Nextflow’s scalability for efficient processing of large BCF files.
- **Data Processing**: Handles BCF format, a standard in genomics, with robust filtering logic.
- **Educational Value**: Teaches students Rust, Nextflow, and parallel processing in a bioinformatics context, aligning with genomic analysis workflows.
- **Scalability**: Suitable for HPC environments, where large genomic datasets require distributed processing.

### Next Steps / Improvements
- **Dynamic Chunking**: Read chunks from a file (`chunks_list.txt`) instead of hardcoding them in `main.nf`.
- **Additional Filters**: Extend `bcf_filter_tool` to support more filtering criteria (e.g., allele frequency, genotype quality).
- **Metrics Reporting**: Add variant statistics (e.g., number of filtered variants per chunk) to output files.
- **Real Data Integration**: Use a larger, real BCF dataset (e.g., from 1000 Genomes) to test scalability.
- **Visualization**: Plot filtering statistics (e.g., QUAL distribution) using Python/Matplotlib.
- **Pipeline Enhancements**: Add error recovery and resume functionality in Nextflow.
- **Output Formats**: Support VCF output alongside BCF for broader compatibility.
