# 5.5. Variant Calling and Genotyping
## Experiment 5.5.1: Parallel Variant Calling for Genomic Hypotheses
This project provides a Rust program to perform parallel variant calling on genomic data, processing variant hypotheses against pileup data to compute likelihoods of alternative alleles. The program is designed for educational purposes, demonstrating parallel processing and chunked data handling in bioinformatics using Rust’s performance capabilities. It processes pileup data (base and quality information) and variant hypotheses from JSON files, computing naive log-likelihoods for each hypothesis. The implementation leverages parallel computing with Rayon and chunked processing to ensure scalability, making it suitable for analyzing large-scale biological datasets in high-performance computing (HPC) environments, although adapted for local execution in this experiment.
## Breakdown of the Project
### Dependencies
The project consists of a single Rust program (`naive_variant_caller`). Dependencies are:
**Rust Program** (`naive_variant_caller/Cargo.toml`):
```
[dependencies]
anyhow = "1.0"
rayon = "1.7"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
clap = { version = "4.0", features = ["derive"] }
```

- **anyhow**: Simplifies error handling for robust execution.
- **rayon**: Enables parallel processing of variant hypotheses within chunks.
- **serde**: Provides serialization/deserialization for JSON input/output.
- **serde**_json: Handles JSON file parsing and writing.
- **clap**: Parses command-line arguments for input files and parameters.

**Feature Flags**:The program has no feature flags, using a single workflow with configurable parameters:

- Default: ``chunk_size=10_000``, `output_dir="partial_variants"`, `merged_output="merged_variants.json"`.
- Data: Expects `pileup.json` (base/quality data) and `hypotheses.json` (variant hypotheses), producing partial and merged JSON files.

The program has one main process:

- naive_variant_caller: Processes hypotheses in chunks, computes likelihoods, writes partial results, and merges them into a final JSON file.

## Implementation Modules
### Rust Program (`naive_variant_caller/src/main.rs`):
### Data Structure:

- Defines `BaseQuality` for base and quality (`{base: char, quality: u8}`).
- Defines `VariantHypothesis` for variant hypotheses (`{chrom: String, position: u64, ref_base: char, alt_base: char}`).
- Defines `VariantSite` for results (`{chrom: String, position: u64, ref_base: char, alt_base: char, likelihood: f64}`).
- Defines `PartialResult` for partial results (`{sites: Vec<VariantSite>}`).
- Uses `HashMap<(String, u64), Vec<BaseQuality>>` for pileup data, with keys parsed from strings like `"chr1:100"`.

### Likelihood Calculation:

- Implements `naive_likelihood` to compute naive `log-likelihoods` using natural log probabilities.
- Interprets base quality as error probability (e.g., `10^(-quality/10))`, computing match probability (`1.0 - p_error`) for reference or alternative base, and error probability (`p_error`) for others.

### Chunk Processing:

- Processes all chunks sequentially in `main`, dividing hypotheses into chunks of `chunk_size`.
- Processes hypotheses in parallel with Rayon’s `par_iter()` in `compute_variants`, ensuring thread-safe data access.
- Writes partial results to `partial_variants_chunk_<index>.json`.

### Merging:

- Implements `merge_variant_sites` to combine partial JSON files into `merged_variants.json`.
- Dynamically reads chunk files from `output_dir using read_dir`.

### Main Function:

- Parses command-line arguments with `clap`.
- Loads pileup.json and `hypotheses.json`.
- Processes all chunks and merges results in a single run.
- Outputs JSON files totaling ~1-2 KB for small test datasets (e.g., 4 hypotheses across 2 chunks).

## Program Execution
The program processes pileup and hypothesis data to produce partial and merged variant call files in JSON format. Steps to execute:
**Prerequisites**:

Install Rust/Cargo:
`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env`

Setup: Ensure `naive_variant_caller/` (containing `Cargo.toml, src/main.rs, pileup.json, hypotheses.json, hpc_job.slurm`) is in the project directory.

Run Program:
`cd ~/BGVR/chapter_05/t1/naive_variant_caller
cargo build --release
bash hpc_job.slurm`

Example Output 
Chunk 0 processed (2 hypotheses). Partial results saved at "partial_variants/partial_variants_chunk_0.json"
Successfully merged 4 variants into "merged_variants.json"

Generated Files:

partial_variants/partial_variants_chunk_0.json: Partial results for chunk 0 (~500 bytes).
partial_variants/partial_variants_chunk_1.json: Partial results for chunk 1 (~500 bytes).
merged_variants.json: Merged variant calls (~1 KB).
variant_caller_*.out, variant_caller_*.err: Log files from hpc_job.slurm.
Located in ~/BGVR/chapter_05/t1/naive_variant_caller/.

Docker Alternative:
docker build -t naive_variant_caller .
mkdir -p partial_variants
docker run -v $(pwd)/partial_variants:/app/partial_variants naive_variant_caller

Note:The hpc_job.slurm script specifies a --chunk-size of 2, but the provided main.rs does not support the --chunk-index argument expected by earlier experiment versions. The output suggests only one chunk was processed, indicating a possible mismatch or partial execution. To align with the expected output (2 chunks), ensure input files contain at least 4 hypotheses.
Why is this Project Important?

Bioinformatics Workflow: Demonstrates variant calling, critical for identifying genomic variants in sequencing studies.
Parallel Processing: Leverages Rayon for efficient hypothesis processing, showcasing Rust’s concurrency.
Data Processing: Introduces chunked processing to manage memory for large datasets, relevant to HPC pipelines.
Educational Value: Teaches students Rust, parallel computing, and JSON-based data handling in a bioinformatics context, aligning with Chapter 5’s focus on variant analysis.
Scalability: Produces output suitable for integration into larger HPC pipelines, where partial results can be merged or analyzed further.

Next Steps / Improvements

Input Flexibility: Support BAM or VCF inputs for pileup and hypothesis data.
Parameterization: Add command-line flags for likelihood model parameters or output formats.
Pipeline Integration: Wrap in Nextflow or Snakemake for automated HPC workflows.
Likelihood Models: Replace naive_likelihood with Bayesian or machine learning-based models.
Performance: Test with larger datasets (e.g., millions of hypotheses) to evaluate scaling.
Output Format: Add VCF output for compatibility with standard genomics tools.
Chunk Indexing: Add --chunk-index support to align with earlier experiment versions, enabling SLURM array jobs.


