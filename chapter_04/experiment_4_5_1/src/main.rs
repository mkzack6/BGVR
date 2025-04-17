use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

/// Represents an exon or genomic segment in a splicing graph.
/// Fields are placeholders for future boundary validation.
#[derive(Debug, Clone, PartialEq)]
struct ExonSegment {
    _start: usize, // Planned: exon start position
    _end: usize,   // Planned: exon end position
}

/// Represents alignment data for one read.
#[derive(Debug, Clone)]
struct Alignment {
    _chrom: String,
    start: usize,
    end: usize,
}

/// A graph representing exons and splice junctions.
/// `adjacency` maps exon coordinates (start, end) to adjacent exons and their coverage.
#[derive(Debug, Default, Clone)]
struct SplicingGraph {
    adjacency: HashMap<(usize, usize), Vec<(ExonSegment, u64)>>,
}

impl SplicingGraph {
    /// Adds or updates a junction with dynamic coverage.
    fn add_junction(&mut self, exon_key: (usize, usize), target_exon: ExonSegment) {
        let entry = self.adjacency.entry(exon_key).or_default();
        if let Some((_, cov)) = entry.iter_mut().find(|(exon, _)| exon._start == target_exon._start && exon._end == target_exon._end) {
            *cov += 1; // Increment coverage
        } else {
            entry.push((target_exon, 1)); // New junction
        }
    }

    /// Merges another SplicingGraph into the current one.
    fn merge(&mut self, other: SplicingGraph) {
        for (key, edges) in other.adjacency {
            self.adjacency.entry(key).or_default().extend(edges);
        }
    }
}

/// Processes a slice of `Alignment` entries to produce a local splicing graph.
fn process_alignment_chunk(batch: &[Alignment]) -> SplicingGraph {
    let mut local_graph = SplicingGraph::default();
    for align in batch {
        if align.start >= align.end {
            continue; // Skip invalid alignments
        }
        let exon_key = (align.start, align.end);
        let target_exon = ExonSegment {
            _start: align.start, // TODO: Derive from CIGAR or annotations
            _end: align.end,
        };
        local_graph.add_junction(exon_key, target_exon);
    }
    local_graph
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let all_alignments = vec![
        Alignment {
            _chrom: "chr1".to_string(),
            start: 100,
            end: 200,
        },
        Alignment {
            _chrom: "chr1".to_string(),
            start: 150,
            end: 300,
        },
        Alignment {
            _chrom: "chr2".to_string(),
            start: 500,
            end: 700,
        },
    ];

    if all_alignments.is_empty() {
        return Err("No alignments provided".into());
    }

    // Chunk size balances memory and parallelism
    let chunk_size = (all_alignments.len() / rayon::current_num_threads()).max(1);
    let chunks: Vec<_> = all_alignments.chunks(chunk_size).map(|c| c.to_vec()).collect();

    let final_graph = chunks
        .into_par_iter()
        .map(|batch| process_alignment_chunk(&batch))
        .reduce(
            SplicingGraph::default,
            |mut acc, local_graph| {
                acc.merge(local_graph);
                acc
            },
        );

    let out_file = File::create("partial_splicing_graph.bin")?;
    let mut writer = BufWriter::new(out_file);
    writer.write_all(b"Splicing graph\n")?;
    writer.write_all(format!("{:#?}", final_graph).as_bytes())?;

    println!("Splicing graph has been successfully written.");

    Ok(())
}
