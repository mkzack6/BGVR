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

-Data Structure:

--Defines `ChromCoverage` to store chromosome names and coverage data as `Vec<f64>`.

--Hardcodes sample data (e.g., `chr2` with 8 positions).

-Smoothing:

--Implements `smooth_coverage` using a prefix-sum approach for O(n) rolling-window averaging.

--Reduces noise in coverage profiles when `do_smooth=true`.

-Peak Calling:

--Uses `local_peak_call` to slide a window, computing mean coverage and marking peaks above `threshold`.

--Processes chromosomes in parallel with `rayon`’s `into_par_iter()`, ensuring thread-safe data handling.

-Output:

--Writes peaks to `partial_peaks.bed` in a BED-like format (`chrom\tpos\tvalue`).

--Uses `BufWriter` with `flush()` to ensure complete file output.

-Main Function:

--Initializes data, calls `call_peaks_and_smooth`, and saves results.

Outputs ~104 bytes for `chr2` peaks (e.g., positions 0–7, values ~4.458–6.611).
### Program Execution
The program processes hardcoded genomic coverage data to produce a peak file in a simplified BED-like format. Steps to execute:

**Prerequisites:**

Install Rust/Cargo:
`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env`
