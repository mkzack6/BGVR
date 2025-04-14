#!/usr/bin/env nextflow

/*
  A minimal Nextflow pipeline that:
    1) Compiles the Rust program (process 'compile')
    2) Executes the Rust binary on the provided FASTQ data (process 'analysis')

  Usage:
    nextflow run main.nf --fastq /path/to/reads.fastq
*/

nextflow.enable.dsl=2

params.fastq = 'example.fastq'
params.outdir = 'results'

process compile {
    tag "Compiling code"
    executor 'local'

    input:
    path project_dir

    output:
    path "${project_dir}/target/release/debruijn_bloom"

    """
    echo "Compiling the Rust code..."
    cd ${project_dir}
    cargo build --release
    """
}

process analysis {
    tag "DeBruijn & Bloom"
    executor 'local'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path debruijn_bloom
    path fastq_file

    output:
    path 'graph.json'
    path 'bloom.json'

    """
    echo "Running the Rust analysis on ${fastq_file}..."
    ./debruijn_bloom --fastq ${fastq_file} --kmer 5 --outdir ${params.outdir}
    cp ${params.outdir}/graph.json ./
    cp ${params.outdir}/bloom.json ./
    """
}

workflow {
    project_ch = Channel.fromPath('debruijn_bloom', type: 'dir', checkIfExists: true)
    fastq_ch = Channel.fromPath(params.fastq, checkIfExists: true)
    compile_out = compile(project_ch)
    analysis(compile_out, fastq_ch)
}