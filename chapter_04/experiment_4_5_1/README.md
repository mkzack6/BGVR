# 4.5. Transcriptomics and Alternative Splicing Algorithms
## Experiment 4.5.1: Parallel Splicing Graph Construction

This project provides a Rust program to construct splicing graphs from RNA-seq alignment data, representing gene transcripts by connecting exons and splice junctions. The program is designed for educational purposes, demonstrating parallel computing in bioinformatics and Rust’s performance capabilities. It processes alignment data in parallel, builds partial splicing graphs, merges them, and outputs the result as a text file. The implementation emphasizes computational efficiency and scalability, suitable for large-scale transcriptomics datasets.
## Breakdown of the Project
### Dependencies
The project consists of a single Rust program (`parallel_splicing_graph`). Dependencies are defined in Cargo.toml:

Rust Program (`parallel_splicing_graph/Cargo.toml`):
`rayon = "1.10.0"`: Enables parallel processing of alignment chunks.
Standard library (`std::collections, std::fs, std::io`): Provides hash maps, file handling, and buffered writing.

### Feature Flags
The Rust program has no feature flags, using a single workflow to process alignments with a configurable chunk size:

Default: `chunk_size` dynamically set based on alignment count and CPU threads, outputs to `partial_splicing_graph.bin`.
Customizable: Accepts a custom `chunk_size` via code modification (not command-line).

### Implementation Modules

Rust Program (`parallel_splicing_graph/src/main.rs`):
**Data Structures:**
`ExonSegment`: Represents an exon with start and end coordinates (placeholders for future validation).
`Alignment`: Stores read alignment data (chromosome, start, end).
`SplicingGraph`: A hash map mapping exon coordinates to adjacent exons and coverage.

**Chunk Processing:**
Splits alignments into chunks using dynamic sizing (`len / num_threads`).
Processes chunks in parallel with `rayon::into_par_iter` to build local splicing graphs.

**Graph Merging**:
Merges partial graphs using a `reduce` operation, combining adjacency maps.

**Output**:
Writes the graph as a debug-formatted text file (`partial_splicing_graph.bin`).
Uses `BufWriter` with `flush` and `sync_all` for reliable file writes.

**Main Function**:
Validates input, processes alignments, merges graphs, and saves the output.

## Program Execution
The program processes alignment data to produce a splicing graph, outputting a text file. Steps to execute:
Prerequisites

`Install Rust/Cargo:curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env`

### Setup

Ensure `Cargo.toml` and `src/main.rs` are in the project directory (`parallel_splicing_graph/`).

### Run Program
`cd ~/BGVR/chapter_04/experiment_4_2_2
cargo run --release`

### Example Output
Splicing graph has been successfully written.

### Generated Files:

`partial_splicing_graph.bin`: Text file containing the debug-formatted splicing graph (~1KB for 3 alignments).
Located in the project root.

## Why is this Project Important?

* **Bioinformatics Workflow**: Demonstrates splicing graph construction, critical for transcriptomics and isoform detection.
* **Parallel Processing**: Leverages `rayon` for efficient graph building, showcasing Rust’s concurrency.
* **Data Structures**: Uses hash maps for graph representation, relevant to genomic data.
* **Educational Value**: Teaches students Rust, parallel computing, and bioinformatics in a reproducible pipeline.

## Next Steps / Improvements

* **Real Data**: Process actual RNA-seq alignments (e.g., BAM files) using CIGAR strings.
* **Output Formats**: Support standard formats like GFA or JSON for interoperability.
* **Validation**: Add checks for exon overlaps and chromosome consistency.
* **Optimization**: Use sparse data structures for large graphs to reduce memory usage.
* **Benchmarking**: Profile performance for larger datasets (e.g., millions of alignments).
* **Integration**: Link with transcript assemblers to validate splicing patterns.


