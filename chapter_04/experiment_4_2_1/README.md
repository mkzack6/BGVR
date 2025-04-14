# 4.2. Graph-Based Models for Gene Regulatory Networks (GRNs)
## Experiment 4.2.1: Parallel Correlation Matrix Pipeline

This pipeline uses Nextflow to run a Rust binary (`parallel_correlation`) that generates a correlation adjacency matrix for a specified number of genes and samples.

## Prerequisites
- Nextflow (`version 24.10.5` or later)
- Rust (for compiling `parallel_correlation`)
- Linux environment (tested on Ubuntu)

## Directory Structure
```
experiment_4_2_1/
├── parallel_correlation/
│   ├── src/
│   │   └── main.rs
│   └── Cargo.toml
├── pipeline.nf
└── partial_adjacency.bin (generated output)
```
