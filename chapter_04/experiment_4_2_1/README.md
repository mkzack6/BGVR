# 4.2. Graph-Based Models for Gene Regulatory Networks (GRNs)
## Experiment 4.2.1: Parallel Correlation Matrix Pipeline

This project provides a Nextflow pipeline that compiles and runs a Rust program to generate a correlation adjacency matrix for simulated gene expression data. The pipeline is designed for educational purposes, demonstrating parallel computing in bioinformatics, Rust’s performance capabilities, and workflow orchestration with Nextflow. It generates random gene expression data, computes pairwise correlations in parallel, and outputs the resulting matrix as a binary file. The implementation emphasizes computational efficiency and scalability, suitable for analyzing high-dimensional biological datasets.
## Breakdown of the Project

#### Dependencies
The project consists of a Nextflow script (`pipeline.nf`) and a Rust program (`parallel_correlation`). Dependencies are:

**Nextflow Pipeline (`pipeline.nf`)**:
- **Nextflow**: Orchestrates compilation and execution of the Rust binary.
- **Rust/Cargo**: Required to compile the Rust code.
- No additional Nextflow plugins are needed; runs locally with executor `local`.

**Rust Program (`parallel_correlation/Cargo.toml`)**:
```toml
[dependencies]
rayon = "1.10.0"
ndarray = "0.16.1"
rand = "=0.9.0"
rand_core = "=0.9.0"
rand_chacha = "=0.9.0"
```
**rayon:** Enables parallel computation of pairwise correlations.

**ndarray:** Provides efficient matrix operations for data generation and correlation calculations.

**rand, rand_core, rand_chacha:** Generate high-quality random numbers for synthetic gene expression data.

### Feature Flags

The Rust program has no feature flags, using a single workflow to generate and process data with configurable parameters:

Default: `--num-genes 1000`, `--num-samples 50`, outputs to `partial_adjacency.bin`.

Customizable: Accepts `--num-genes <size>`, `--num-samples <size>`, `--output <file>` via command-line (passed by Nextflow).

The Nextflow pipeline has one process:

**buildAdjacencyMatrix:** Compiles and runs the Rust binary to produce the correlation matrix.
## Implementation Modules
### Nextflow Pipeline (`pipeline.nf`):

-buildAdjacencyMatrix Process:

--Inputs: Number of genes (`val genes`), number of samples (`val samples`), output file name (`val output_file`).

--Outputs: Correlation matrix (`partial_adjacency.bin`).

--Copies the `parallel_correlation/` directory to the work directory.

--Runs `cargo build --release` to compile the Rust binary.

--Executes the binary with `--num-genes`, `--num-samples`, and `--output`` to generate the matrix.

--Includes debugging commands (`echo`, `ls`) to verify file creation.

--Uses `publishDir` to copy the output to the project root and `sync` for reliable file writes..

-Workflow: Executes `buildAdjacencyMatrix` with DSL2 for modularity.

### Rust Program (`parallel_correlation/src/main.rs`):

-Command-Line Parsing:

--Uses `std::env` to parse `--num-genes`, `--num-samples`, and `--output` manually, ensuring lightweight execution.

-Matrix Generation:

--Generates a random `num_genes x num_samples` matrix using `rand::thread_rng()` and `gen_range()` (e.g., `0.0..1000.0` for 1000x50).

--Stores data in `ndarray::Array2<f64>` for efficient access.

-Correlation Calculation:

--Computes pairwise Pearson correlations in parallel using `rayon`, targeting the upper triangle and mirroring results.

--Uses `Arc<Mutex<Array2<f64>>>` for thread-safe updates to the `num_genes x num_genes` matrix.

-Output:

--Writes the matrix as a binary file (`partial_adjacency.bin`) using `f64` values.

--Includes `writer.flush()` and `file.sync_all()` to ensure file integrity.

--Verifies file existence with `Path::new().exists()` to confirm successful creation.

-Main Function:

--Parses arguments, generates data, computes correlations, and saves the matrix.

--Outputs matrix size (~8MB for 1000x1000).
## Program Execution
The pipeline compiles the Rust program and runs it to produce a correlation matrix, outputting a binary file. Steps to execute:

### Prerequisites:

**Install Rust/Cargo:**
`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env`
**Install Nextflow:**
`curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
chmod +x ~/bin/nextflow`
**Setup:**
Ensure `pipeline.nf` and `parallel_correlation/` (containing `Cargo.toml` and `src/main.rs`) are in the project directory.

**Run Pipeline:**
`cd ~/BGVR/chapter_04/experiment_4_2_1
nextflow run pipeline.nf --num_genes 1000 --num_samples 50`
**Example Output:**
```
N E X T F L O W  ~  version 24.10.5
Launching `pipeline.nf` [maniac_heyrovsky] DSL2 - revision: dcd56f9f39
executor >  local (1)
[5b/bdfff5] buildAdjacencyMatrix (Build Adjacency Matrix) [100%] 1 of 1 ✔
```
**Generated Files:**

-`partial_adjacency.bin`: Binary file containing a 1000x1000 `f64` correlation matrix (~8MB).

-Located in `cd ~/BGVR/chapter_04/experiment_4_2_1` (via `publishDir`) and `work/5b/bdfff5.../`.

## Why is this Project Important?
Bioinformatics Workflow: Demonstrates correlation analysis, critical for gene expression studies and co-expression networks.

Pipeline Orchestration: Introduces Nextflow for automating Rust compilation and execution, scalable to larger datasets.

Parallel Processing: Leverages `rayon` for efficient matrix operations, showcasing Rust’s concurrency.

Data Structures: Uses dense matrices for correlation storage, relevant to high-dimensional biological data.

Educational Value: Teaches students Rust, Nextflow, and parallel computing in a reproducible bioinformatics pipeline.
## Next Steps / Improvements
Real Data: Replace random data with actual gene expression datasets (e.g., RNA-seq counts).

Matrix Optimization: Use sparse matrices for large-scale correlations to reduce memory usage.

Output Formats: Add support for text-based outputs (e.g., CSV) or compressed binaries.

Distributed Processing: Extend pipeline with cluster executors (e.g., SLURM) for big datasets.

Integration: Link with `debruijn_bloom` to analyze correlations in k-mer graphs.

Benchmarking: Profile performance for larger inputs (e.g., `--num_genes 10000`).
