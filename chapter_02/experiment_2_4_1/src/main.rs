use std::collections::HashMap;
use rayon::prelude::*;
use bio::io::fasta;

/// Represents a node in the pangenome graph, storing a k-mer label.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct PGNode {
    kmer: String,
}

/// An adjacency list mapping each node to its successor nodes.
type PGGraph = HashMap<PGNode, Vec<PGNode>>;

/// Builds a pangenome graph by merging k-mers from multiple FASTA files.
/// Each file corresponds to one haplotype/genome. Overlapping k-mers
/// are treated as edges, and equivalent k-mers across files are merged.
fn build_pangenome_graph(k: usize, fasta_files: &[&str]) -> PGGraph {
    // Read each FASTA in parallel, then build a partial k-mer graph for each.
    let partial_graphs: Vec<PGGraph> = fasta_files
        .par_iter()
        .map(|path| {
            let reader = fasta::Reader::from_file(path)
                .unwrap_or_else(|_| panic!("Cannot open FASTA file: {}", path));
            
            let mut local_graph = PGGraph::new();
            // Process each record in the FASTA.
            for record in reader.records() {
                let rec = record.expect("Invalid FASTA record");
                let seq_bytes = rec.seq();
                // If the sequence is shorter than k, skip it.
                if seq_bytes.len() < k {
                    continue;
                }
                let seq_str = String::from_utf8(seq_bytes.to_vec())
                    .expect("Non-UTF8 sequence data");
                
                // Build a local adjacency for k-mers.
                for i in 0..seq_str.len().saturating_sub(k) {
                    let node_kmer = &seq_str[i..i + k];
                    let edge_kmer = &seq_str[i + 1..i + k + 1];

                    let node = PGNode {
                        kmer: node_kmer.to_owned(),
                    };
                    let edge_node = PGNode {
                        kmer: edge_kmer.to_owned(),
                    };

                    local_graph
                        .entry(node)
                        .or_insert_with(Vec::new)
                        .push(edge_node);
                }
            }
            local_graph
        })
        .collect();

    // Reduce step: merge the partial graphs into one global pangenome graph.
    // We use .into_par_iter() so we can supply an identity closure (|| PGGraph::new())
    // and a combining closure |acc, local|.
    partial_graphs
        .into_par_iter()
        .reduce(
            || PGGraph::new(),
            |mut acc, local| {
                for (node, successors) in local {
                    acc.entry(node)
                        .or_insert_with(Vec::new)
                        .extend(successors);
                }
                acc
            },
        )
}

fn main() {
    // Example usage: merging multiple haplotypes, each in its own FASTA file.
    let haplotypes = &["src/haplotype1.fasta", "src/haplotype2.fasta", "src/haplotype3.fasta"];
    let k = 21;
    let pangenome_graph = build_pangenome_graph(k, haplotypes);

    println!(
        "Constructed a pangenome graph with {} nodes.",
        pangenome_graph.len()
    );
    // Print a small subset of the graph
    for (node, edges) in pangenome_graph.iter().take(5) {
        println!("Node: {} -> {:?}", node.kmer, 
                 edges.iter().map(|e| &e.kmer).collect::<Vec<_>>());
    }
}
