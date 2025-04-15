# 4.5. Transcriptomics and Alternative Splicing Algorithms
## Experiment 4.5.1: Parallel Splicing Graph Construction
This project provides a Rust program to build and merge partial splicing graphs, representing gene transcripts by connecting exons and splice junctions from alignment data, such as in RNA-seq workflows. The program is designed for educational purposes, demonstrating parallel computing in bioinformatics and Rust’s performance capabilities. It processes read alignments, constructs local splicing graphs in chunks, merges them in parallel, and outputs the resulting graph as a serialized file. The implementation emphasizes computational efficiency and scalability, suitable for analyzing large-scale transcriptomic datasets.
## Breakdown of the Project
### Dependencies
The project consists of a single Rust program (splicing_graph). Dependencies are:
Rust Program (splicing_graph/Cargo.toml):
[dependencies]
rayon = "1.10.0"

rayon: Enables parallel processing of alignment chunks to construct and merge splicing graphs efficiently.

Feature Flags:

The Rust program has no feature flags, using a single workflow with hardcoded parameters:
Default: chunk_size=100, outputs to partial_splicing_graph.bin.
Data: Hardcoded alignments for chromosomes (e.g., chr1, chr2), producing a graph with exon junctions.

The program has one main process:

splicing_graph: Processes alignments, builds local graphs, merges them, and writes the result.

Implementation Modules
Rust Program (splicing_graph/src/main.rs):

Data Structures:
Defines ExonSegment to store exon boundaries (_start, _end), with fields underscored as unused.
Defines Alignment to hold read alignment data (_chrom, start, end).
Defines SplicingGraph as a HashMap mapping exon keys (start, end) to vectors of target exons and coverage.

Graph Construction:
Implements process_alignment_chunk to create a local SplicingGraph from a batch of alignments.
Adds junctions using add_junction, linking exon keys to target exons with coverage.

Parallel Processing:
Splits alignments into chunks using chunks(chunk_size).
Uses rayon’s into_par_iter() and reduce to build and merge graphs in parallel.

Output:
Writes the graph to partial_splicing_graph.bin as a debug string, prefixed with a header.
Uses BufWriter with flush() to ensure complete file output.

Main Function:
Initializes hardcoded alignments, processes them in chunks, merges graphs, and saves the result.
Outputs ~320 bytes for a graph with three exon junctions (e.g., (100, 200) to _start: 150, _end: 300).

Program Execution
The program processes hardcoded alignment data to produce a serialized splicing graph in a text-based format. Steps to execute:
Prerequisites:

Install Rust/Cargo:
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

Setup: Ensure splicing_graph/ (containing Cargo.toml and src/main.rs) is in the project directory:
~/BGVR/chapter_04/experiment_4_5_1/
├── splicing_graph/
│   ├── Cargo.toml
│   ├── src/
│   │   └── main.rs

Run Program:
cd ~/BGVR/chapter_04/experiment_4_5_1/splicing_graph
cargo build --release
./target/release/splicing_graph

Example Output:
Splicing graph has been successfully written.
-rw-r--r-- 1 zack zack 320 Apr 15 02:00 partial_splicing_graph.bin

Content of partial_splicing_graph.bin:
Serialized splicing graph example
SplicingGraph {
    adjacency: {
        (100, 200): [
            (
                ExonSegment {
                    _start: 150,
                    _end: 300
                },
                1
            )
        ],
        (150, 300): [
            (
                ExonSegment {
                    _start: 200,
                    _end: 400
                },
                1
            )
        ],
        (500, 700): [
            (
                ExonSegment {
                    _start: 550,
                    _end: 800
                },
                1
            )
        ]
    }
}

Generated Files:

partial_splicing_graph.bin: Text file containing a serialized splicing graph (~320 bytes).
Located in ~/BGVR/chapter_04/experiment_4_5_1/splicing_graph/.

Why is this Project Important?

Bioinformatics Workflow: Demonstrates splicing graph construction, critical for inferring transcript structures in RNA-seq studies.
Parallel Processing: Leverages rayon for efficient processing of alignment data, showcasing Rust’s concurrency.
Data Structures: Uses a HashMap-based graph to store exon junctions, relevant to transcriptomic analyses.
Educational Value: Teaches students Rust and parallel computing in a bioinformatics context, aligning with Chapter 4’s focus on high-performance workflows.
Scalability: Produces output suitable for merging in larger HPC pipelines, enabling analysis of extensive sequencing datasets.

Next Steps / Improvements

File Output Reliability: Add writer.flush() to BufWriter to ensure all data is written to disk, preventing potential incomplete output.
Input Flexibility: Support reading alignments from files (e.g., BAM, SAM) using libraries like noodles.
Parameterization: Introduce command-line arguments for chunk_size and output file name.
Pipeline Integration: Wrap in a Nextflow or Snakemake pipeline for automated HPC workflows.
Graph Optimization: Aggregate coverage values during merging to reduce redundancy.
Integration: Link with experiment_4_2_1 (correlation analysis) or experiment_4_4_1 (peak calling) to combine transcriptomic and genomic analyses.


