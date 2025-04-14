use rayon::prelude::*;
use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::BufWriter;
use std::path::PathBuf;
use std::{env, error::Error};

#[derive(Debug, Serialize, Deserialize)]
struct Node {
    kmer: String,
    edges: HashSet<String>,
}

#[derive(Debug, Serialize, Deserialize)]
struct BloomFilter {
    bits: Vec<bool>,
    num_hashes: usize,
    size: usize,
}

impl BloomFilter {
    fn new(size: usize, num_hashes: usize) -> Self {
        BloomFilter {
            bits: vec![false; size],
            num_hashes,
            size,
        }
    }

    fn insert(&self, item: &str) {
        for seed in 0..self.num_hashes {
            let h = self.hash_with_seed(item, seed) % self.size;
            unsafe {
                let bits_ptr = self.bits.as_ptr() as *mut bool;
                *bits_ptr.add(h) = true;
            }
        }
    }

    fn contains(&self, item: &str) -> bool {
        for seed in 0..self.num_hashes {
            let h = self.hash_with_seed(item, seed) % self.size;
            if !self.bits[h] {
                return false;
            }
        }
        true
    }

    fn hash_with_seed(&self, item: &str, seed: usize) -> usize {
        use std::hash::{Hash, Hasher};
        use std::collections::hash_map::DefaultHasher;

        let mut hasher = DefaultHasher::new();
        item.hash(&mut hasher);
        seed.hash(&mut hasher);
        hasher.finish() as usize
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = env::args();
    let _bin = args.next();

    let mut fastq_path = String::from("example.fastq");
    let mut kmer_size = 31;
    let mut outdir = String::from("results");

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--fastq" => fastq_path = args.next().unwrap(),
            "--kmer" => {
                let val = args.next().unwrap();
                kmer_size = val.parse().unwrap_or(31);
            }
            "--outdir" => outdir = args.next().unwrap(),
            _ => {}
        }
    }

    fs::create_dir_all(&outdir)?;

    let mut reader = parse_fastx_file(PathBuf::from(&fastq_path))?;
    println!("Reading from FASTQ: {}, k-mer={}", fastq_path, kmer_size);

    let mut reads: Vec<String> = Vec::new();
    while let Some(record) = reader.next() {
        let seqrec = record?;
        reads.push(String::from_utf8_lossy(seqrec.seq().as_ref()).to_string());
    }
    println!("Loaded {} reads. Building de Bruijn graph & Bloom filter...", reads.len());

    let nodes_map = build_debruijn(&reads, kmer_size);

    let distinct_kmers_count = nodes_map.len();
    let bf_size = (distinct_kmers_count * 10).max(1_000_000);
    let bloom = BloomFilter::new(bf_size, 3);

    nodes_map
        .par_iter()
        .for_each(|(kmer, _node)| {
            bloom.insert(kmer);
        });

    let test_kmer = "ACGT";
    println!("Bloom filter contains '{}'? {}", test_kmer, bloom.contains(test_kmer));

    let graph_path = format!("{}/graph.json", outdir);
    let bloom_path = format!("{}/bloom.json", outdir);

    let graph_file = File::create(&graph_path)?;
    serde_json::to_writer(BufWriter::new(graph_file), &nodes_map)?;

    let bloom_file = File::create(&bloom_path)?;
    serde_json::to_writer(BufWriter::new(bloom_file), &bloom)?;

    println!("De Bruijn graph saved to {}", graph_path);
    println!("Bloom filter saved to {}", bloom_path);
    Ok(())
}

fn build_debruijn(reads: &[String], k: usize) -> HashMap<String, Node> {
    let kmers = reads
        .par_iter()
        .flat_map(|read| {
            let mut local_kmers = Vec::new();
            if read.len() >= k {
                for i in 0..=read.len() - k {
                    local_kmers.push(&read[i..i + k]);
                }
            }
            local_kmers
        })
        .collect::<Vec<&str>>();

    let mut map: HashMap<String, Node> = HashMap::new();
    for kmer_str in kmers.iter() {
        map.entry((*kmer_str).to_string())
            .or_insert(Node {
                kmer: (*kmer_str).to_string(),
                edges: HashSet::new(),
            });
    }

    reads.iter().for_each(|read| {
        if read.len() < k + 1 {
            return;
        }
        for i in 0..=read.len() - (k + 1) {
            let k1 = &read[i..i + k];
            let k2 = &read[i + 1..i + 1 + k];
            if let Some(node) = map.get_mut(k1) {
                node.edges.insert(k2.to_string());
            }
        }
    });

    map
}
