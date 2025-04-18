use anyhow::Result;
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::fs::{File, create_dir_all};
use std::io::BufWriter;
use clap::Parser;

#[derive(Parser)]
struct Args {
    #[arg(long)]
    pileup_input: String,
    #[arg(long)]
    hypotheses_input: String,
    #[arg(long, default_value_t = 1000)]
    chunk_size: usize,
    #[arg(long, default_value_t = 0)]
    chunk_index: usize,
    #[arg(long)]
    output_dir: String,
    #[arg(long)]
    merged_output: String,
}

#[derive(Serialize, Deserialize)]
struct PileupEntry {
    base: String,
    quality: u8,
}

#[derive(Serialize, Deserialize)]
struct Hypothesis {
    chrom: String,
    position: u64,
    ref_base: String,
    alt_base: String,
}

#[derive(Serialize, Deserialize)]
struct VariantSite {
    chrom: String,
    position: u64,
    ref_base: String,
    alt_base: String,
    likelihood: f64,
}

#[derive(Serialize, Deserialize)]
struct ChunkOutput {
    sites: Vec<VariantSite>,
}

fn calculate_likelihood(pileup: &[PileupEntry], ref_base: &str, alt_base: &str) -> f64 {
    let mut log_likelihood = 0.0;
    for entry in pileup {
        let error_prob = 10.0_f64.powf(-(entry.quality as f64) / 10.0);
        let match_prob = 1.0 - error_prob;
        let prob = if entry.base == ref_base {
            match_prob
        } else if entry.base == alt_base {
            error_prob / 3.0
        } else {
            error_prob / 3.0
        };
        log_likelihood += prob.log10();
    }
    log_likelihood
}

fn process_chunk(hypotheses: &[Hypothesis], pileup: &HashMap<String, Vec<PileupEntry>>, chunk_size: usize, chunk_index: usize, output_dir: &str) -> Result<()> {
    let start = chunk_index * chunk_size;
    let end = (start + chunk_size).min(hypotheses.len());
    if start >= hypotheses.len() {
        return Ok(());
    }

    let chunk_hypotheses = &hypotheses[start..end];
    let sites: Vec<VariantSite> = chunk_hypotheses.par_iter().map(|hyp| {
        let key = format!("{}:{}", hyp.chrom, hyp.position);
        static EMPTY: Vec<PileupEntry> = Vec::new();
        let pileup_data = pileup.get(&key).unwrap_or(&EMPTY);
        let likelihood = calculate_likelihood(pileup_data, &hyp.ref_base, &hyp.alt_base);
        VariantSite {
            chrom: hyp.chrom.clone(),
            position: hyp.position,
            ref_base: hyp.ref_base.clone(),
            alt_base: hyp.alt_base.clone(),
            likelihood,
        }
    }).collect();

    let chunk_output = ChunkOutput { sites };
    let chunk_file = format!("{}/partial_variants_chunk_{}.json", output_dir, chunk_index);
    let file = BufWriter::new(File::create(&chunk_file)?);
    serde_json::to_writer(file, &chunk_output)?;
    println!("Chunk {} processed ({} hypotheses). Partial results saved at {:?}", chunk_index, chunk_hypotheses.len(), chunk_file);
    Ok(())
}

fn merge_variants(output_dir: &str, merged_output: &str, num_chunks: usize) -> Result<()> {
    let mut merged_variants = vec![];
    for i in 0..num_chunks {
        let chunk_file = format!("{}/partial_variants_chunk_{}.json", output_dir, i);
        let file = File::open(&chunk_file)?;
        let chunk_output: ChunkOutput = serde_json::from_reader(file)?;
        merged_variants.extend(chunk_output.sites);
    }
    let file = BufWriter::new(File::create(merged_output)?);
    serde_json::to_writer(file, &merged_variants)?;
    println!("Successfully merged {} variants into {:?}", merged_variants.len(), merged_output);
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    create_dir_all(&args.output_dir)?;

    let pileup_file = File::open(&args.pileup_input)?;
    let pileup: HashMap<String, Vec<PileupEntry>> = serde_json::from_reader(pileup_file)?;

    let hypotheses_file = File::open(&args.hypotheses_input)?;
    let hypotheses: Vec<Hypothesis> = serde_json::from_reader(hypotheses_file)?;

    process_chunk(&hypotheses, &pileup, args.chunk_size, args.chunk_index, &args.output_dir)?;

    let num_chunks = (hypotheses.len() + args.chunk_size - 1) / args.chunk_size;
    if args.chunk_index == 0 {
        merge_variants(&args.output_dir, &args.merged_output, num_chunks)?;
    }

    Ok(())
}
