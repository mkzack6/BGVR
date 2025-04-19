# 5.6. Structural Variant Detection
## Experiment 5.6.1: Parallel Breakpoint Detection from Alignment Segments
This project provides a Rust program to detect breakpoints indicative of structural variants (SVs) from alignment segment data, such as those derived from split-read alignments in genomic sequencing. The program is designed for educational purposes, demonstrating parallel processing in bioinformatics using Rust’s high-performance capabilities. It processes alignment segments in JSON format, groups them by read ID, and identifies breakpoints based on consecutive segment properties (e.g., intra-chromosomal gaps or inter-chromosomal jumps). The implementation leverages parallel computing with `rayon` to ensure scalability, making it suitable for analyzing large-scale genomic datasets in high-performance computing (HPC) environments.
## Breakdown of the Project
### Dependencies
The project consists of a single Rust program (`experiment_5_6_1`). Dependencies are defined in `experiment_5_6_1/Cargo.toml`:
[dependencies]
anyhow = "1.0"
rayon = "1.7"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
clap = { version = "4.3", features = ["derive"] }


- **anyhow**: Provides flexible error handling for robust program execution.
- **rayon**: Enables parallel processing of read groups to detect breakpoints efficiently.
- **serde and serde_json**: Support JSON parsing and serialization of alignment segments and breakpoints.
- **clap**: Facilitates command-line argument parsing for input file, chunk size, and output paths.

### Feature Flags
The program has no feature flags, using a single workflow with configurable parameters:

- **Default**: `chunk_size=5000`, `partial_output_dir="partial_breakpoints"`, `merged_output="merged_breakpoints.json"`.
- **Data**: Expects a JSON file (`alignments.json`) containing alignment segments, producing partial and merged breakpoint files in JSON format.

The program has one main process:

- **experiment_5_6_1**: Groups alignment segments by read ID, detects breakpoints, writes partial results per chunk, and merges them into a final output.

### Implementation Modules
Rust Program (`experiment_5_6_1/src/main.rs`):
### Data Structures:

- `AlignmentSegment`: Stores read ID, chromosome, start position, CIGAR string, and orientation (e.g., `{"read_id": "read1", "chrom": "chr1", "start": 1000, "cigar": "50M", "orientation": "+"}`).
- `Breakpoint`: Represents detected breakpoints with read ID, chromosome, position, and SV type (e.g., `"intra-chr"` or `"translocation"`).
- `PartialBreakpoints`: Holds a vector of breakpoints for each chunk, serialized to JSON.

### Breakpoint Detection:

- Implements `detect_breakpoints` to analyze consecutive segments of the same read:
  - For segments on the same chromosome, calculates breakpoint position as `start + CIGAR length` (e.g., `1000 + 50 = 1050` for `50M`).
  - For segments on different chromosomes, marks a translocation at the first segment’s start position.


- Processes read groups in parallel using `rayon`’s `par_iter()`, ensuring thread-safe data handling.

### Chunk Processing:

- Reads segments in chunks (default size 5000) to manage memory for large inputs.
- Groups segments by read ID using a `HashMap`, sorts by start position, and detects breakpoints.

### Input/Output:

- Reads `alignments.json`, supporting both array-based (e.g., `[{...}, {...}]`) and map-based (e.g., `{"segments": [...]}`) JSON formats.
- Writes partial breakpoint files to `partial_breakpoints/partial_breakpoints_{index}.json`.
- Merges partial files into `merged_breakpoints.json` in JSON array format.
- Uses `BufWriter` for efficient file output with proper flushing.

### Main Function:

- Parses command-line arguments with `clap`.
- Loads and parses JSON input, processes chunks, detects breakpoints, and merges results.
- Outputs ~2 breakpoints for the example input (4 segments, 2 read groups), producing files of ~200–300 bytes.

## Program Execution
The program processes alignment segments from `alignments.json` to detect breakpoints and produce JSON output files. Steps to execute:
### Prerequisites

Install Rust/Cargo:`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env`



### Setup

Ensure the `experiment_5_6_1/` directory (containing `Cargo.toml` and `src/main.rs`) is in the project directory:`~/BGVR/chapter_05/experiment_5_6_1/`


Verify `alignments.json` exists in the same directory, e.g.:
```
[
    {"read_id": "read1", "chrom": "chr1", "start": 1000, "cigar": "50M", "orientation": "+"},
    {"read_id": "read1", "chrom": "chr1", "start": 5000, "cigar": "50M", "orientation": "+"},
    {"read_id": "read2", "chrom": "chr2", "start": 2000, "cigar": "50M", "orientation": "-"},
    {"read_id": "read2", "chrom": "chr3", "start": 3000, "cigar": "50M", "orientation": "+"}
]
```


### Run Program
```
cd ~/BGVR/chapter_05/experiment_5_6_1/
cargo build --release
./target/release/experiment_5_6_1 --alignment-input alignments.json
```

### Example Output
`Parsed JSON: Array [Object {"chrom": String("chr1"), "cigar": String("50M"), ...}, ...]
Processed chunk 0 (2 read groups). Wrote partial results to "partial_breakpoints/partial_breakpoints_0.json"
Merged 2 total breakpoints into "merged_breakpoints.json".`

### Generated Files

- **partial_breakpoints/partial_breakpoints_0.json**: JSON file containing breakpoints for the chunk (~200 bytes).
```
{"breakpoints":[
    {"read_id":"read1","chrom":"chr1","pos":1050,"sv_type":"intra-chr"},
    {"read_id":"read2","chrom":"chr2","pos":2000,"sv_type":"translocation"}
]}
```
- **merged_breakpoints.json**: Merged breakpoints in JSON array format (~200 bytes).
```
[
    {"read_id":"read1","chrom":"chr1","pos":1050,"sv_type":"intra-chr"},
    {"read_id":"read2","chrom":"chr2","pos":2000,"sv_type":"translocation"}
]
```

- **Location**: Files are generated in `~/BGVR/chapter_05/experiment_5_6_1/t1/ and ~/BGVR/chapter_05/experiment_5_6_1/partial_breakpoints/`.

## Why is this Project Important?

- **Bioinformatics Workflow**: Demonstrates breakpoint detection for structural variant analysis, critical for understanding genomic rearrangements in sequencing studies.
- **Parallel Processing**: Leverages rayon for efficient processing of read groups, showcasing Rust’s concurrency capabilities.
- **Data Processing**: Handles JSON-based alignment data, a common format in bioinformatics pipelines, with robust error handling.
- **Educational Value**: Teaches students Rust, parallel computing, and JSON parsing in a bioinformatics context, aligning with Chapter 5’s focus on genomic variant detection.
- **Scalability**: Produces output suitable for integration into larger HPC pipelines, where breakpoint files can be further analyzed or visualized.

## Next Steps / Improvements

- **Input Flexibility**: Add support for reading alignment data from BAM or SAM files using crates like rust-htslib.
- **Parameterization**: Introduce command-line arguments for CIGAR parsing rules or SV type definitions.
- **Pipeline Integration**: Wrap in a Nextflow or Snakemake pipeline for automated HPC workflows.
- **Breakpoint Refinement**: Implement filtering of low-confidence breakpoints or merging of nearby breakpoints.
- **Performance**: Test with larger datasets (e.g., thousands of reads) to evaluate parallel scaling and optimize memory usage with streaming deserialization.
- **Output Formats**: Support additional output formats like VCF for compatibility with standard SV analysis tools.
