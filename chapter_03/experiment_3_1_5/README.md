# 3.1 Introduction to Data Structures and Algorithms
**Experiment 3.1.5: Genomic Indexer**

This Rust program demonstrates a distributed genomic indexing workflow using MPI (Message Passing Interface). It simulates indexing a genome across multiple processes (ranks), where each rank builds a partial index, sends it to rank 0, which merges the indexes and saves the result to `output.txt`. The program is designed for educational purposes, illustrating parallel computing concepts in bioinformatics, such as splitting computational tasks across nodes in a high-performance computing (HPC) cluster.

## Purpose

In bioinformatics, indexing a genome (e.g., human DNA) enables fast searches for patterns, like finding genes or aligning sequencing reads. This program mimics that process in a simplified way:
- Each rank (process) indexes a portion of the genome, creating a `PartialIndex`.
- Indexes are serialized and sent to rank 0 using point-to-point communication.
- Rank 0 merges the partial indexes into a global index and saves it to `output.txt`.

The workflow reflects real-world tools like BWA, which index genomes for efficient read alignment.

## Prerequisites

- **Rust**: Install via `rustup` (`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`).
- **OpenMPI**: Install for MPI support (`sudo apt-get install openmpi-bin libopenmpi-dev` on Ubuntu/WSL).
- **Cargo**: Comes with Rust for building the project.
- **Dependencies**: Specified in `Cargo.toml`.

## Project Structure

- `src/main.rs`: The main program implementing the distributed indexing workflow.
- `Cargo.toml`: Defines dependencies and project settings.
- `output.txt`: Generated file containing the merged global index.

## Code Overview

The program uses Rust with the `mpi` crate for parallelism and `bincode`/`serde` for serialization. Key components:

- **`PartialIndex`**: A struct representing a rank’s index chunk, with a rank ID and a `HashMap` simulating DNA pattern positions.
- **`build_local_index`**: Creates a partial index for a rank (e.g., `{"key_rank_1": 1}`).
- **`merge_indexes`**: Combines partial indexes into a global `HashMap` at rank 0.
- **`main`**:
  - Initializes MPI to assign ranks (0–3 for `-np 4`).
  - Each rank builds and serializes its index.
  - Ranks 1–3 send serialized data to rank 0; rank 0 receives and adds its own.
  - Rank 0 deserializes, merges, and saves the global index to `output.txt`.

Here’s the complete code (`src/main.rs`):

```rust
use mpi::traits::*;
use std::collections::HashMap;
use serde::{Serialize, Deserialize};
use bincode;
use std::fs::File;
use std::io::Write;
use std::env;

/// A chunk of a genomic index built by one MPI rank (e.g., a node in an HPC cluster).
/// It’s like a piece of a DNA search engine for part of the genome.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct PartialIndex {
    rank_id: i32,                     // Which rank built this chunk
    index_data: HashMap<String, usize>, // Mini-database of index data
}

/// Builds a partial index for this rank’s slice of the genome.
/// In practice, this might index DNA sequences (e.g., for a suffix array).
fn build_local_index(rank: i32) -> PartialIndex {
    let mut data = HashMap::new();
    data.insert(format!("key_rank_{}", rank), rank as usize); // Simulate index entry
    PartialIndex {
        rank_id: rank,
        index_data: data,
    }
}

/// Merges partial indexes into a global index.
/// Like combining puzzle pieces into a full genome index.
fn merge_indexes(chunks: Vec<PartialIndex>) -> HashMap<String, usize> {
    let mut global_map = HashMap::new();
    for chunk in chunks {
        for (k, v) in chunk.index_data {
            global_map.insert(k, v); // Copy entries to global index
        }
    }
    global_map
}

fn main() {
    // Log MPI environment to confirm communicator setup
    if let Ok(rank) = env::var("OMPI_COMM_WORLD_RANK") {
        println!("OMPI_COMM_WORLD_RANK: {}", rank);
    }
    if let Ok(size) = env::var("OMPI_COMM_WORLD_SIZE") {
        println!("OMPI_COMM_WORLD_SIZE: {}", size);
    }

    // Initialize MPI for parallel processing
    let universe = mpi::initialize().expect("Failed to initialize MPI");
    let world = universe.world();
    let rank = world.rank(); // Rank ID (0, 1, 2, 3 for -np 4)
    let size = world.size(); // Number of ranks (4 for -np 4)
    println!("Rank {} of {} starting", rank, size);

    // Build and serialize this rank’s index
    let local_chunk = build_local_index(rank);
    let serialized_chunk = bincode::serialize(&local_chunk)
        .expect("Failed to serialize partial index");

    // Sync ranks before communication
    world.barrier();

    // Rank 0 collects all indexes; others send
    let mut gathered_chunks = Vec::new();
    if rank == 0 {
        gathered_chunks.push(serialized_chunk); // Rank 0’s own chunk
        for source_rank in 1..size {
            // Receive serialized chunk from each rank
            let (received_chunk, _status) = world
                .process_at_rank(source_rank)
                .receive_vec::<u8>();
            gathered_chunks.push(received_chunk);
        }
    } else {
        // Non-root ranks send to rank 0
        world.process_at_rank(0).send(&serialized_chunk[..]);
    }

    // Sync after communication
    world.barrier();

    // Only rank 0 merges and saves
    if rank == 0 {
        // Deserialize chunks
        let deserialized_chunks: Vec<PartialIndex> = gathered_chunks
            .into_iter()
            .map(|bytes| bincode::deserialize(&bytes).expect("Failed to deserialize"))
            .collect();

        // Merge into global index
        let global_index = merge_indexes(deserialized_chunks);

        // Save to output.txt
        let mut file = File::create("output.txt")
            .expect("Failed to create output.txt");
        writeln!(file, "Final merged index: {:?}", global_index)
            .expect("Failed to write to output.txt");

        // Print for verification
        println!("Final merged index on rank 0: {:?}", global_index);
    }
}
