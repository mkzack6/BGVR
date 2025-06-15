#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Genomic Sparse Index Pipeline
 * 
 * This pipeline processes genomic data files using a Rust-based sparse indexing tool
 * and performs parallel range queries for allele frequency analysis.
 */

// Default parameters
params.input = "genotypes.txt"
params.output_dir = "results"
params.max_index = 10000000
params.num_tasks = 8
params.help = false

// Help message
if (params.help) {
    log.info """
    Genomic Sparse Index Pipeline
    =============================
    
    Usage:
        nextflow run main.nf [options]
    
    Options:
        --input         Input genotype file (default: genotypes.txt)
        --output_dir    Output directory (default: results)
        --max_index     Maximum genomic position index (default: 10000000)
        --num_tasks     Number of parallel query tasks (default: 8)
        --help          Show this help message
    
    Example:
        nextflow run main.nf --input my_genotypes.txt --output_dir my_results --num_tasks 16
    """
    exit 1
}

// Log pipeline parameters
log.info """
Genomic Sparse Index Pipeline
=============================
Input file       : ${params.input}
Output directory : ${params.output_dir}
Max index        : ${params.max_index}
Number of tasks  : ${params.num_tasks}
"""

/*
 * Validate input file
 */
process validateInput {
    tag "validation"
    
    input:
    path input_file
    
    output:
    path input_file, emit: validated_file
    path "validation_report.txt", emit: report
    
    script:
    """
    echo "Validating input file: ${input_file}" > validation_report.txt
    
    if [ ! -f "${input_file}" ]; then
        echo "ERROR: Input file does not exist" >> validation_report.txt
        exit 1
    fi
    
    # Check file format and content
    line_count=\$(wc -l < "${input_file}")
    echo "Total lines: \$line_count" >> validation_report.txt
    
    # Check for required columns in first few non-comment lines
    valid_lines=\$(grep -v '^#' "${input_file}" | head -10 | awk 'NF >= 3 {print \$0}' | wc -l)
    echo "Valid data lines (sample): \$valid_lines" >> validation_report.txt
    
    if [ \$valid_lines -eq 0 ]; then
        echo "ERROR: No valid data lines found (expected format: chr pos count [freq])" >> validation_report.txt
        exit 1
    fi
    
    echo "Validation completed successfully" >> validation_report.txt
    """
}

/*
 * Build sparse index using Rust tool
 */
process buildSparseIndex {
    tag "indexing"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path input_file
    val max_index
    val num_tasks
    
    output:
    path "index_results.txt", emit: results
    path "rust_sparse_index.log", emit: log
    
    script:
    """
    # Ensure the Rust binary is available
    if [ ! -f "./rust_sparse_index" ]; then
        echo "ERROR: rust_sparse_index binary not found in current directory"
        echo "Please compile the Rust code first with: cargo build --release"
        echo "Then copy target/release/rust_sparse_index to the workflow directory"
        exit 1
    fi
    
    # Run the sparse index builder with logging
    ./rust_sparse_index \\
        --input "${input_file}" \\
        --output "index_results.txt" \\
        --max-index ${max_index} \\
        --num-tasks ${num_tasks} \\
        --verbose > rust_sparse_index.log 2>&1
    
    # Check if the process completed successfully
    if [ \$? -ne 0 ]; then
        echo "ERROR: Sparse index building failed"
        cat rust_sparse_index.log
        exit 1
    fi
    
    echo "Sparse index building completed successfully" >> rust_sparse_index.log
    """
}

/*
 * Generate summary statistics
 */
process generateSummary {
    tag "summary"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path results_file
    path validation_report
    path process_log
    
    output:
    path "pipeline_summary.html", emit: summary
    path "pipeline_stats.txt", emit: stats
    
    script:
    """
    # Generate text summary
    echo "Genomic Sparse Index Pipeline Summary" > pipeline_stats.txt
    echo "====================================" >> pipeline_stats.txt
    echo "Run date: \$(date)" >> pipeline_stats.txt
    echo "Input parameters:" >> pipeline_stats.txt
    echo "  - Input file: ${params.input}" >> pipeline_stats.txt
    echo "  - Max index: ${params.max_index}" >> pipeline_stats.txt
    echo "  - Query tasks: ${params.num_tasks}" >> pipeline_stats.txt
    echo "" >> pipeline_stats.txt
    
    # Add validation results
    echo "Validation Results:" >> pipeline_stats.txt
    cat "${validation_report}" >> pipeline_stats.txt
    echo "" >> pipeline_stats.txt
    
    # Add processing results
    echo "Processing Results:" >> pipeline_stats.txt
    if [ -f "${results_file}" ]; then
        result_lines=\$(wc -l < "${results_file}")
        echo "  - Results file lines: \$result_lines" >> pipeline_stats.txt
        echo "  - Query results:" >> pipeline_stats.txt
        grep -v '^#' "${results_file}" | head -5 >> pipeline_stats.txt
        if [ \$(grep -v '^#' "${results_file}" | wc -l) -gt 5 ]; then
            echo "  ... (showing first 5 results)" >> pipeline_stats.txt
        fi
    else
        echo "  - ERROR: Results file not found" >> pipeline_stats.txt
    fi
    
    # Generate HTML summary
    cat > pipeline_summary.html << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>Genomic Sparse Index Pipeline Summary</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { background-color: #f0f0f0; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; }
        .code { background-color: #f8f8f8; padding: 10px; border-left: 3px solid #007acc; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <div class="header">
        <h1>Genomic Sparse Index Pipeline</h1>
        <p><strong>Run Date:</strong> \$(date)</p>
        <p><strong>Status:</strong> Completed Successfully</p>
    </div>
    
    <div class="section">
        <h2>Pipeline Parameters</h2>
        <table>
            <tr><th>Parameter</th><th>Value</th></tr>
            <tr><td>Input File</td><td>${params.input}</td></tr>
            <tr><td>Output Directory</td><td>${params.output_dir}</td></tr>
            <tr><td>Max Index</td><td>${params.max_index}</td></tr>
            <tr><td>Query Tasks</td><td>${params.num_tasks}</td></tr>
        </table>
    </div>
    
    <div class="section">
        <h2>Results Files Generated</h2>
        <ul>
            <li><strong>index_results.txt</strong> - Main results with range query summaries</li>
            <li><strong>pipeline_stats.txt</strong> - Detailed statistics and metrics</li>
            <li><strong>rust_sparse_index.log</strong> - Processing logs</li>
            <li><strong>pipeline_summary.html</strong> - This summary report</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>Next Steps</h2>
        <p>The sparse genomic index has been successfully built. You can now:</p>
        <ul>
            <li>Examine the range query results in <code>index_results.txt</code></li>
            <li>Use the Rust binary directly for custom queries</li>
            <li>Integrate the results into downstream genomic analysis workflows</li>
        </ul>
    </div>
</body>
</html>
EOF
    """
}

/*
 * Main workflow
 */
workflow {
    // Check if input file exists
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Validate input
    validateInput(input_ch)
    
    // Build sparse index
    buildSparseIndex(
        validateInput.out.validated_file,
        params.max_index,
        params.num_tasks
    )
    
    // Generate summary
    generateSummary(
        buildSparseIndex.out.results,
        validateInput.out.report,
        buildSparseIndex.out.log
    )
    
    // Final output
    buildSparseIndex.out.results.view { "Results: $it" }
    generateSummary.out.summary.view { "Summary: $it" }
}
