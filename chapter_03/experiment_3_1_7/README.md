# 3.1. Introduction to Data Structures and Algorithms 
## Experiment 3.1.7: De Bruijn Graph and Bloom Filter Pipeline

This project provides a Nextflow pipeline that compiles and runs a Rust program to process FASTQ files, constructing a de Bruijn graph and a Bloom filter from k-mers. The pipeline is designed for educational purposes, demonstrating bioinformatics workflows, parallel processing with Rust, and pipeline orchestration with Nextflow. It reads FASTQ sequences, extracts k-mers, builds a graph of k-mer adjacencies, stores k-mers in a Bloom filter for membership queries, and outputs results as JSON files. The implementation emphasizes memory efficiency and scalability, suitable for genomic data analysis.

## Breakdown of the Project

### 1. Dependencies
The project consists of a Nextflow script (`main.nf`) and a Rust program (`debruijn_bloom`). Dependencies are:

- **Nextflow Pipeline (`main.nf`)**:
  - **Nextflow**: Orchestrates compilation and analysis processes.
  - **Rust/Cargo**: Required for compiling the Rust code.
  - No additional Nextflow plugins are needed; runs locally with `executor 'local'`.

- **Rust Program (`debruijn_bloom/Cargo.toml`)**:
  ```toml 
  [dependencies]
  rayon = "1.7"
  needletail = "0.6"
  serde = { version = "1", features = ["derive"] }
  serde_json = "1"
  ```
rayon: Enables parallel k-mer extraction and Bloom filter insertion.

needletail: Parses FASTQ files efficiently with streaming.

serde, serde_json: Serializes graph and Bloom filter to JSON.

### Feature Flags
The Rust program has no feature flags, using a single workflow to process FASTQ files with configurable parameters:

Default: Processes `example.fastq` with `kmer=5`, outputs to `results/`.

Customizable: Accepts `--fastq <path>`, `--kmer <size>`, `--outdir <dir>` via command-line (passed by Nextflow).

The Nextflow pipeline has two processes:

compile: Builds the Rust binary.

analysis: Runs the binary on the FASTQ file.

## Implementation Modules
Nextflow Pipeline (`main.nf`):

compile Process:

Inputs: `debruijn_bloom` project directory.

Outputs: Compiled binary (`target/release/debruijn_bloom`).

Runs `cargo build --release` in the staged project directory.

analysis Process:

Inputs: Compiled binary, FASTQ file.

Outputs: `graph.json` (de Bruijn graph), `bloom.json` (Bloom filter).

Executes the binary with `--kmer 5` for short reads.

Workflow: Chains compilation to analysis, using DSL2 for modularity.

Rust Program (`debruijn_bloom/src/main.rs`):

Node Struct:

Represents a k-mer with edges to adjacent k-mers.

Serialized to JSON: `{"kmer":"ACGTA","edges":["CGTAC", ...]}`.

BloomFilter Struct:

Stores k-mers in a bit vector with multiple hash functions.

Supports parallel insertion using `rayon`.

Serialized to JSON: `{"bits":[false,...],"num_hashes":3,"size":1000000}`.

build_debruijn Function:

Extracts k-mers in parallel (O(n) per read, n=read length).

Builds a `HashMap<String, Node>` for the graph.

main Function:

Reads FASTQ with `needletail`.

Constructs graph and Bloom filter.

Writes outputs to `results/graph.json` and `results/bloom.json`

## Program Execution
The pipeline compiles the Rust program and runs it on a FASTQ file, producing JSON outputs. Steps to execute:

### Prerequisites:
Install Rust/Cargo:
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```
Install Nextflow:
```
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
chmod +x ~/bin/nextflow
```
### Setup
Ensure `main.nf`, `debruijn_bloom/`, and `example.fastq` are in the project directory.
`debruijn_bloom/` contains `Cargo.toml` and `src/main.rs`.
### Run Pipeline:
```
cd ~/(replace with project_directory)
nextflow run main.nf --fastq example.fastq
```
### Example Output
```
N E X T F L O W  ~  version 24.10.5
Launching `main.nf` [id] DSL2 - revision: ...
executor >  local (2)
[Compiling code] Compiling the Rust code...
[DeBruijn & Bloom] Running the Rust analysis on example.fastq...
Reading from FASTQ: example.fastq, k-mer=5
Loaded 2 reads. Building de Bruijn graph & Bloom filter...
Bloom filter contains 'ACGT'? true
De Bruijn graph saved to results/graph.json
Bloom filter saved to results/bloom.json
```
Files: `results/graph.json`, `results/bloom.json`

## Why is this Project Important?
Bioinformatics Workflow: Demonstrates k-mer analysis, a foundation for genome assembly and sequence alignment.

Pipeline Orchestration: Introduces Nextflow for automating compilation and analysis, scalable to HPC clusters.

Parallel Processing: Uses rayon for efficient k-mer extraction, showcasing Rust’s concurrency.

Data Structures: Combines de Bruijn graphs (graph theory) with Bloom filters (probabilistic membership), key to genomics.

Educational Value: Teaches students Rust, Nextflow, and genomic data processing in a reproducible pipeline.
## Next Steps / Improvements
Larger Datasets: Test with real FASTQ files (e.g., Illumina reads).

K-mer Optimization: Add variable k-mer sizes or canonical k-mers.

Bloom Filter: Replace unsafe code with thread-safe bitvec crate.

Distributed Processing: Extend pipeline with MPI or cluster executors (e.g., SLURM).

Integration: Combine with fasta_parser to preprocess FASTA inputs.

Benchmarking: Profile performance on large genomic datasets.
