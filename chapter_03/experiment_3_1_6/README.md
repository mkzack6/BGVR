# 3.1: Introduction to Data Structures and Algorithms
## Experiment 3.1.6: FASTA Parser

This Rust project demonstrates efficient parsing and analysis of FASTA files using the `needletail` crate. It provides a streamlined interface for reading nucleotide sequences, computing statistics like GC content, and outputting results in JSON format. The program is designed for educational purposes, illustrating memory-efficient bioinformatics workflows suitable for large genomic datasets.

## Breakdown of the Project

### 1. Dependencies (Cargo.toml)
The project uses:

- **needletail**: For fast, memory-efficient parsing of FASTA/FASTQ files.
- **serde**: For serializing/deserializing data structures.
- **serde_json**: For outputting results in JSON format.

```toml
[dependencies]
needletail = "0.5.1"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
```
### Feature Flags
This project does not use feature flags, as it targets a single sequential parsing workflow. All processing occurs in a unified code path, with optional command-line input for the FASTA file path.

Default: Parses a FASTA file (defaults to example.fasta), computes GC content, and writes JSON output.
### Implementation Modules
The implementation is contained in a single module (main.rs) with key components:

-SeqRecord Struct:

  --Stores a sequence’s ID, length, and GC content.

  --Uses serde for JSON serialization.

  --Example: {"id":"seq1","seq_len":9,"gc_content":0.4444444444444444}.

-calc_gc_content Function:

  --Computes the fraction of G and C nucleotides in a sequence.

  --Operates on a borrowed Cow<'_, [u8]> slice from needletail for efficiency.

  --Time complexity: O(ℓ) per sequence, where ℓ is sequence length.

-main Function:

  --Reads a FASTA file using needletail::parse_fastx_file.

  --Streams records to minimize memory usage.

  --Collects SeqRecords and writes to output.json.
### Program Execution
The `main()` function processes a FASTA file, computes statistics, and outputs JSON. It accepts a file path via command-line argument or defaults to `example.fasta`.

Build and run commands:
```
# Build the project
cargo build --release

# Run with default file (example.fasta)
./target/release/fasta_parser

# Run with custom file
./target/release/fasta_parser custom.fasta
```
Example Output:

`Processed 3 sequences`

## Why is this Project Important?
Memory Efficiency: Uses `needletail`’s streaming to handle large genomic files without loading them fully into memory.

Bioinformatics Relevance: Introduces sequence analysis, a key step in read alignment and variant calling.

Rust Benefits: Ensures safety and performance with zero-cost abstractions and `Cow` slices.

Scalability: Linear O(n × ℓ) complexity suits large datasets, with potential for parallelization.

Educational Value: Teaches students file parsing, data serialization, and genomic statistics.

## Next Steps / Improvements
Parallel Processing: Add `rayon` to parallelize GC content calculations for multicore CPUs.

Extended Statistics: Compute additional metrics (e.g., AT content, sequence complexity).

File Format Support: Extend to FASTQ or compressed files.

Integration: Combine with tools like `genomic_indexer` for full pipelines.

Benchmarking: Measure performance on real genomic datasets.
