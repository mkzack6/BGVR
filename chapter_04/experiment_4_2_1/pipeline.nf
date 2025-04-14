#!/usr/bin/env nextflow

/*
 * Usage example:
 *   nextflow run pipeline.nf --num_genes 1000 --num_samples 50
 */

nextflow.enable.dsl=2

params.num_genes = 1000
params.num_samples = 50
params.output_file = "partial_adjacency.bin"

process buildAdjacencyMatrix {
    tag "Build Adjacency Matrix"
    executor 'local'
    stageOutMode 'copy'
    publishDir "${baseDir}", mode: 'copy', pattern: "${params.output_file}"

    input:
    val genes
    val samples
    val output_file

    output:
    path "${output_file}"

    """
    echo "Work directory: \$PWD"
    echo "Listing work directory contents before execution"
    ls -l \$PWD
    cp -r /home/zack/BGVR/chapter_04/experiment_4_2_1/parallel_correlation .
    cd parallel_correlation
    cargo build --release
    echo "Running binary in \$PWD"
    ./target/release/parallel_correlation \
        --num-genes ${genes} \
        --num-samples ${samples} \
        --output \$PWD/../${output_file}
    sync
    echo "Checking output in work directory"
    ls -l \$PWD/../${output_file} || echo "Output file not found in work directory"
    echo "Listing final work directory contents"
    ls -l \$PWD/..
    """
}

workflow {
    buildAdjacencyMatrix(params.num_genes, params.num_samples, params.output_file)
}