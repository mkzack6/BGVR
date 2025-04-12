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


