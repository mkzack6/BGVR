# 4.4. Epigenomic Data Integration and Algorithms
## Experiment 4.4.1: Parallel Peak Calling for Genomic Coverage Data

This project provides a Rust program to identify peaks (regions of high signal intensity) in genomic coverage data, such as from ChIP-seq or ATAC-seq experiments. The program is designed for educational purposes, demonstrating parallel processing in bioinformatics using Rust’s performance capabilities. It processes chromosomal coverage data, optionally smooths it with a rolling average, and identifies peaks where the mean coverage exceeds a threshold. The implementation leverages parallel computing to ensure scalability, making it suitable for analyzing large-scale biological datasets in high-performance computing (HPC) environments.

## Breakdown of the Project

### Dependencies
The project consists of a single Rust program (`peak_calling`). Dependencies are:

**Rust Program (`peak_calling/Cargo.toml`)**:
```toml
[dependencies]
rayon = "1.10.0"
```
**rayon:** Enables parallel processing of multiple chromosomes to identify peaks efficiently.
### Feature Flags:

The program has no feature flags, using a single workflow with hardcoded parameters:

Default: `window=3`, `threshold=3.0`, `do_smooth=true`, outputs to `partial_peaks.bed`.

Data: Hardcoded coverage for chromosomes (e.g., `chr2`), producing peaks in a BED-like format.

The program has one main process:

**peak_calling:** Smooths coverage data, identifies peaks, and writes results to a file.
### Implementation Modules
Rust Program (`peak_calling/src/main.rs`):

**Data Structure**:
- Defines `ChromCoverage` to store chromosome names and coverage data as `Vec<f64>`.
- Hardcodes sample data (e.g., `chr2` with 8 positions).
**Smoothing**:
- Implements `smooth_coverage` using a prefix-sum approach for O(n) rolling-window averaging.
- Reduces noise in coverage profiles when `do_smooth=true`.
**Peak Calling**:
- Uses `local_peak_call` to slide a window, computing mean coverage and marking peaks above `threshold`.
- Processes chromosomes in parallel with `rayon`’s `into_par_iter()`, ensuring thread-safe data handling.
**Output**:
- Writes peaks to `partial_peaks.bed` in a BED-like format (`chrom\tpos\tvalue`).
- Uses `BufWriter` with `flush()` to ensure complete file output.
**Main Function**:
- Initializes data, calls `call_peaks_and_smooth`, and saves results.
Outputs ~104 bytes for `chr2` peaks (e.g., positions 0–7, values ~4.458–6.611).
### Program Execution
The program processes hardcoded genomic coverage data to produce a peak file in a simplified BED-like format. Steps to execute:

**Prerequisites:**

Install Rust/Cargo:
`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env`
**Setup:** Ensure `peak_calling/` (containing `Cargo.toml` and `src/main.rs`) is in the project directory:
**Run Program:**
`cd ~/BGVR/chapter_04/experiment_4_4_1/peak_calling
cargo build --release
./target/release/peak_calling`
**Example Output:**
`-rw-r--r-- 1 zack zack 104 Apr 15 01:41 partial_peaks.bed`
**Generated Files:**

`partial_peaks.bed`: Tab-separated file containing chromosome, position, and peak value (~104 bytes).

Located in `~/BGVR/chapter_04/experiment_4_4_1/peak_calling/`.

## Why is this Project Important?

**Bioinformatics Workflow:** Demonstrates peak calling, critical for identifying regulatory regions in genomic studies like ChIP-seq.

**Parallel Processing:** Leverages `rayon` for efficient analysis of chromosomal data, showcasing Rust’s concurrency.

**Data Processing:** Introduces smoothing techniques to handle noisy biological signals, relevant to real-world datasets.

**Educational Value:** Teaches students Rust and parallel computing in a bioinformatics context, aligning with Chapter 4’s focus on high-performance workflows.

**Scalability:** Produces output suitable for integration into larger HPC pipelines, where peak files can be merged or analyzed further.

## Next Steps / Improvements

**Input Flexibility:** Add support for reading coverage data from files (e.g., CSV, BAM, bigWig).

**Parameterization:** Introduce command-line arguments for `window`, `threshold`, and `do_smooth`.

**Pipeline Integration:** Wrap in a Nextflow or Snakemake pipeline for automated HPC workflows.

**Peak Refinement:** Implement merging of nearby peaks or filtering of spurious calls.

**Performance:** Test with larger datasets to evaluate parallel scaling (e.g., multiple chromosomes).
