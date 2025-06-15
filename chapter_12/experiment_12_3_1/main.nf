#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Enhanced Genomic Algorithms Pipeline
 * 
 * This pipeline performs comprehensive genomic analysis including:
 * - Allele frequency calculation
 * - Linkage disequilibrium analysis
 * - HMM-based haplotype phasing
 * - Quality control metrics
 */

// Default parameters
params.input = "genotypes.csv"
params.output_dir = "genomic_results"
params.format = "txt"
params.hmm_states = 3
params.seed = 42
params.compute_ld = false
params.max_ld_pairs = 1000
params.phase_all = false
params.min_maf = 0.01
params.help = false

// Help message
if (params.help) {
    log.info """
    Enhanced Genomic Algorithms Pipeline
    ===================================
    
    This pipeline performs comprehensive genomic analysis using advanced algorithms
    implemented in Rust for high performance and memory efficiency.
    
    Usage:
        nextflow run main.nf [options]
    
    Required inputs:
        --input         Input genotype CSV file (required)
    
    Analysis options:
        --output_dir    Output directory (default: genomic_results)
        --format        Output format: txt, json, csv (default: txt)
        --hmm_states    Number of HMM states for phasing (default: 3)
        --seed          Random seed for reproducibility (default: 42)
        --compute_ld    Enable linkage disequilibrium computation (default: false)
        --max_ld_pairs  Maximum LD pairs to compute (default: 1000)
        --phase_all     Phase all individuals vs sample (default: false)
        --min_maf       Minimum minor allele frequency (default: 0.01)
        --help          Show this help message
    
    Examples:
        # Basic analysis
        nextflow run main.nf --input genotypes.csv
        
        # Comprehensive analysis with LD and phasing
        nextflow run main.nf \\
            --input large_cohort.csv \\
            --output_dir full_analysis \\
            --compute_ld \\
            --phase_all \\
            --max_ld_pairs 5000 \\
            --format json
        
        # Population genetics study
        nextflow run main.nf \\
            --input population_data.csv \\
            --min_maf 0.05 \\
            --hmm_states 5 \\
            --compute_ld \\
            --format csv
    """
    exit 1
}

// Log pipeline parameters
log.info """
Enhanced Genomic Algorithms Pipeline
===================================
Input file       : ${params.input}
Output directory : ${params.output_dir}
Output format    : ${params.format}
HMM states       : ${params.hmm_states}
Random seed      : ${params.seed}
Compute LD       : ${params.compute_ld}
Max LD pairs     : ${params.max_ld_pairs}
Phase all        : ${params.phase_all}
Min MAF          : ${params.min_maf}
"""

/*
 * Validate and prepare input data
 */
process validateGenomicInput {
    tag "validation"
    
    input:
    path input_file
    
    output:
    path input_file, emit: validated_file
    path "validation_report.txt", emit: report
    path "data_summary.json", emit: summary
    
    script:
    """
    echo "Validating genomic input file: ${input_file}" > validation_report.txt
    echo "Validation started: \$(date)" >> validation_report.txt
    
    if [ ! -f "${input_file}" ]; then
        echo "ERROR: Input file does not exist" >> validation_report.txt
        exit 1
    fi
    
    # Check file format
    if [[ "${input_file}" != *.csv ]]; then
        echo "WARNING: Expected CSV file, got ${input_file}" >> validation_report.txt
    fi
    
    # Count lines and check format
    total_lines=\$(wc -l < "${input_file}")
    echo "Total lines: \$total_lines" >> validation_report.txt
    
    # Check for valid CSV format
    first_data_line=\$(grep -v '^#' "${input_file}" | head -1)
    if [ -n "\$first_data_line" ]; then
        num_fields=\$(echo "\$first_data_line" | tr ',' '\\n' | wc -l)
        echo "Number of variants (fields): \$num_fields" >> validation_report.txt
        
        # Validate genotype values
        invalid_genotypes=\$(grep -v '^#' "${input_file}" | head -100 | grep -o '[^,]' | grep -v '[0-2]' | wc -l)
        if [ \$invalid_genotypes -gt 0 ]; then
            echo "WARNING: Found \$invalid_genotypes invalid genotype values in first 100 lines" >> validation_report.txt
        fi
    else
        echo "ERROR: No valid data lines found" >> validation_report.txt
        exit 1
    fi
    
    # Generate data summary JSON
    cat > data_summary.json << EOF
{
    "input_file": "${input_file}",
    "total_lines": \$total_lines,
    "estimated_variants": \$num_fields,
    "validation_date": "\$(date -Iseconds)",
    "status": "valid"
}
EOF
    
    echo "Validation completed successfully" >> validation_report.txt
    echo "Data appears to be in valid CSV genotype format" >> validation_report.txt
    """
}

/*
 * Run core genomic algorithms
 */
process runGenomicAlgorithms {
    tag "analysis"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path input_file
    val format
    val hmm_states
    val seed
    val compute_ld
    val max_ld_pairs
    val phase_all
    val min_maf
    
    output:
    path "genomic_results.${format}", emit: results
    path "analysis.log", emit: log
    path "performance_metrics.txt", emit: metrics
    
    script:
    def ld_flag = compute_ld ? "--compute-ld" : ""
    def phase_flag = phase_all ? "--phase-all" : ""
    
    """
    # Ensure the Rust binary is available
    if [ ! -f "./rust_core_algorithms" ]; then
        echo "ERROR: rust_core_algorithms binary not found"
        echo "Please compile with: cargo build --release"
        echo "Then copy target/release/rust_core_algorithms to workflow directory"
        exit 1
    fi
    
    echo "Starting genomic analysis: \$(date)" > analysis.log
    echo "Parameters:" >> analysis.log
    echo "  Input: ${input_file}" >> analysis.log
    echo "  Format: ${format}" >> analysis.log
    echo "  HMM states: ${hmm_states}" >> analysis.log
    echo "  Seed: ${seed}" >> analysis.log
    echo "  Compute LD: ${compute_ld}" >> analysis.log
    echo "  Max LD pairs: ${max_ld_pairs}" >> analysis.log
    echo "  Phase all: ${phase_all}" >> analysis.log
    echo "  Min MAF: ${min_maf}" >> analysis.log
    echo "" >> analysis.log
    
    # Record start time and memory
    start_time=\$(date +%s)
    echo "Start time: \$(date)" > performance_metrics.txt
    echo "Initial memory usage:" >> performance_metrics.txt
    free -h >> performance_metrics.txt
    echo "" >> performance_metrics.txt
    
    # Run the analysis with resource monitoring
    /usr/bin/time -v ./rust_core_algorithms \\
        --input "${input_file}" \\
        --output "genomic_results.${format}" \\
        --format "${format}" \\
        --hmm-states ${hmm_states} \\
        --seed ${seed} \\
        ${ld_flag} \\
        --max-ld-pairs ${max_ld_pairs} \\
        ${phase_flag} \\
        --min-maf ${min_maf} \\
        --verbose >> analysis.log 2>&1
    
    exit_code=\$?
    end_time=\$(date +%s)
    duration=\$((end_time - start_time))
    
    # Record performance metrics
    echo "End time: \$(date)" >> performance_metrics.txt
    echo "Total duration: \${duration} seconds" >> performance_metrics.txt
    echo "Final memory usage:" >> performance_metrics.txt
    free -h >> performance_metrics.txt
    echo "" >> performance_metrics.txt
    echo "Exit code: \$exit_code" >> performance_metrics.txt
    
    if [ \$exit_code -ne 0 ]; then
        echo "ERROR: Analysis failed with exit code \$exit_code" >> analysis.log
        cat analysis.log
        exit \$exit_code
    fi
    
    echo "Analysis completed successfully in \${duration} seconds" >> analysis.log
    
    # Validate output file was created
    if [ ! -f "genomic_results.${format}" ]; then
        echo "ERROR: Output file not created" >> analysis.log
        exit 1
    fi
    
    output_size=\$(stat -c%s "genomic_results.${format}")
    echo "Output file size: \$output_size bytes" >> performance_metrics.txt
    
    echo "Genomic analysis pipeline completed successfully" >> analysis.log
    """
}

/*
 * Generate comprehensive analysis report
 */
process generateAnalysisReport {
    tag "reporting"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path results_file
    path analysis_log
    path performance_metrics
    path validation_report
    path data_summary
    
    output:
    path "analysis_report.html", emit: html_report
    path "summary_statistics.txt", emit: summary
    path "pipeline_provenance.json", emit: provenance
    
    script:
    """
    # Generate summary statistics
    echo "Enhanced Genomic Analysis Pipeline - Summary Report" > summary_statistics.txt
    echo "=================================================" >> summary_statistics.txt
    echo "Run date: \$(date)" >> summary_statistics.txt
    echo "Pipeline version: 0.2.0" >> summary_statistics.txt
    echo "" >> summary_statistics.txt
    
    echo "Input Parameters:" >> summary_statistics.txt
    echo "  Input file: ${params.input}" >> summary_statistics.txt
    echo "  Output format: ${params.format}" >> summary_statistics.txt
    echo "  HMM states: ${params.hmm_states}" >> summary_statistics.txt
    echo "  Random seed: ${params.seed}" >> summary_statistics.txt
    echo "  Compute LD: ${params.compute_ld}" >> summary_statistics.txt
    echo "  Max LD pairs: ${params.max_ld_pairs}" >> summary_statistics.txt
    echo "  Phase all: ${params.phase_all}" >> summary_statistics.txt
    echo "  Min MAF: ${params.min_maf}" >> summary_statistics.txt
    echo "" >> summary_statistics.txt
    
    # Add validation results
    echo "Data Validation:" >> summary_statistics.txt
    cat "${validation_report}" | grep -E "(Total lines|Number of variants|WARNING|ERROR)" >> summary_statistics.txt
    echo "" >> summary_statistics.txt
    
    # Add performance metrics
    echo "Performance Metrics:" >> summary_statistics.txt
    cat "${performance_metrics}" | grep -E "(duration|memory|size)" >> summary_statistics.txt
    echo "" >> summary_statistics.txt
    
    # Extract key results if possible
    if [ -f "${results_file}" ]; then
        echo "Analysis Results:" >> summary_statistics.txt
        if [[ "${results_file}" == *.txt ]]; then
            grep -E "(Total|Mean|Variants)" "${results_file}" | head -10 >> summary_statistics.txt
        else
            echo "  Results file: ${results_file}" >> summary_statistics.txt
            echo "  File size: \$(stat -c%s "${results_file}") bytes" >> summary_statistics.txt
        fi
    fi
    
    # Generate HTML report
    cat > analysis_report.html << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>Enhanced Genomic Analysis Report</title>
    <style>
        body { 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            margin: 40px; 
            background-color: #f8f9fa;
        }
        .container { 
            max-width: 1200px; 
            margin: 0 auto; 
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .header { 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px; 
            border-radius: 10px; 
            margin-bottom: 30px;
        }
        .header h1 { margin: 0; font-size: 2.5em; }
        .header p { margin: 10px 0 0 0; font-size: 1.2em; opacity: 0.9; }
        .section { 
            margin: 30px 0; 
            padding: 20px;
            border-left: 4px solid #667eea;
            background-color: #f8f9fa;
        }
        .section h2 { 
            color: #2c3e50; 
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 10px;
        }
        .grid { 
            display: grid; 
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); 
            gap: 20px; 
            margin: 20px 0;
        }
        .card { 
            background: white; 
            padding: 20px; 
            border-radius: 8px; 
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            border: 1px solid #e9ecef;
        }
        .card h3 { 
            color: #495057; 
            margin-top: 0; 
            border-bottom: 1px solid #dee2e6;
            padding-bottom: 10px;
        }
        .metric { 
            display: flex; 
            justify-content: space-between; 
            margin: 10px 0;
            padding: 8px 0;
            border-bottom: 1px dotted #dee2e6;
        }
        .metric:last-child { border-bottom: none; }
        .metric-value { 
            font-weight: bold; 
            color: #28a745; 
        }
        .code { 
            background-color: #f1f3f4; 
            padding: 15px; 
            border-radius: 5px; 
            font-family: 'Consolas', 'Monaco', monospace;
            overflow-x: auto;
            border: 1px solid #e1e4e8;
        }
        .status-success { color: #28a745; font-weight: bold; }
        .status-warning { color: #ffc107; font-weight: bold; }
        .status-error { color: #dc3545; font-weight: bold; }
        .footer {
            margin-top: 40px;
            padding: 20px;
            background-color: #f8f9fa;
            border-radius: 5px;
            text-align: center;
            color: #6c757d;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Enhanced Genomic Analysis</h1>
            <p>Comprehensive population genetics pipeline results</p>
            <p><strong>Generated:</strong> \$(date)</p>
        </div>
        
        <div class="section">
            <h2>📊 Analysis Overview</h2>
            <div class="grid">
                <div class="card">
                    <h3>Input Parameters</h3>
                    <div class="metric">
                        <span>Input File:</span>
                        <span class="metric-value">${params.input}</span>
                    </div>
                    <div class="metric">
                        <span>Output Format:</span>
                        <span class="metric-value">${params.format}</span>
                    </div>
                    <div class="metric">
                        <span>HMM States:</span>
                        <span class="metric-value">${params.hmm_states}</span>
                    </div>
                    <div class="metric">
                        <span>Random Seed:</span>
                        <span class="metric-value">${params.seed}</span>
                    </div>
                </div>
                
                <div class="card">
                    <h3>Analysis Options</h3>
                    <div class="metric">
                        <span>Compute LD:</span>
                        <span class="metric-value">${params.compute_ld ? 'Yes' : 'No'}</span>
                    </div>
                    <div class="metric">
                        <span>Max LD Pairs:</span>
                        <span class="metric-value">${params.max_ld_pairs}</span>
                    </div>
                    <div class="metric">
                        <span>Phase All:</span>
                        <span class="metric-value">${params.phase_all ? 'Yes' : 'No'}</span>
                    </div>
                    <div class="metric">
                        <span>Min MAF:</span>
                        <span class="metric-value">${params.min_maf}</span>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="section">
            <h2>📁 Output Files</h2>
            <div class="card">
                <h3>Generated Files</h3>
                <ul>
                    <li><strong>genomic_results.${params.format}</strong> - Main analysis results</li>
                    <li><strong>analysis.log</strong> - Detailed processing logs</li>
                    <li><strong>performance_metrics.txt</strong> - Performance and resource usage</li>
                    <li><strong>summary_statistics.txt</strong> - Key summary statistics</li>
                    <li><strong>analysis_report.html</strong> - This comprehensive report</li>
                </ul>
            </div>
        </div>
        
        <div class="section">
            <h2>🔬 Analysis Components</h2>
            <div class="grid">
                <div class="card">
                    <h3>Allele Frequency Analysis</h3>
                    <p>Computed population allele frequencies for all variants with quality metrics including call rates and minor allele frequencies.</p>
                </div>
                
                <div class="card">
                    <h3>Linkage Disequilibrium</h3>
                    <p>${params.compute_ld ? 'Calculated LD metrics (r² and D\') between variant pairs to identify population structure and recombination patterns.' : 'LD analysis was not performed (use --compute_ld to enable).'}</p>
                </div>
                
                <div class="card">
                    <h3>HMM Haplotype Phasing</h3>
                    <p>Applied Hidden Markov Model to infer haplotype phase information using ${params.hmm_states} hidden states.</p>
                </div>
                
                <div class="card">
                    <h3>Quality Control</h3>
                    <p>Applied stringent QC filters including MAF ≥ ${params.min_maf} and call rate thresholds.</p>
                </div>
            </div>
        </div>
        
        <div class="section">
            <h2>🚀 Next Steps</h2>
            <div class="card">
                <h3>Recommended Follow-up Analysis</h3>
                <ul>
                    <li><strong>Population Structure:</strong> Use results for PCA or ADMIXTURE analysis</li>
                    <li><strong>Association Studies:</strong> Incorporate phased haplotypes in GWAS</li>
                    <li><strong>Selection Analysis:</strong> Examine LD patterns for signatures of selection</li>
                    <li><strong>Imputation:</strong> Use phasing results for genotype imputation</li>
                </ul>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>Enhanced Genomic Algorithms Pipeline v0.2.0</strong></p>
            <p>Powered by Rust + Nextflow | High-performance population genomics</p>
        </div>
    </div>
</body>
</html>
EOF
    
    # Generate provenance information
    cat > pipeline_provenance.json << EOF
{
    "pipeline_name": "Enhanced Genomic Algorithms",
    "pipeline_version": "0.2.0",
    "execution_date": "\$(date -Iseconds)",
    "parameters": {
        "input": "${params.input}",
        "output_dir": "${params.output_dir}",
        "format": "${params.format}",
        "hmm_states": ${params.hmm_states},
        "seed": ${params.seed},
        "compute_ld": ${params.compute_ld},
        "max_ld_pairs": ${params.max_ld_pairs},
        "phase_all": ${params.phase_all},
        "min_maf": ${params.min_maf}
    },
    "workflow_engine": "Nextflow",
    "compute_environment": "\$(uname -a)",
    "output_files": [
        "genomic_results.${params.format}",
        "analysis.log",
        "performance_metrics.txt",
        "summary_statistics.txt",
        "analysis_report.html"
    ]
}
EOF
    """
}

/*
 * Main workflow
 */
workflow {
    // Validate input file exists
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Validate input data
    validateGenomicInput(input_ch)
    
    // Run core genomic algorithms
    runGenomicAlgorithms(
        validateGenomicInput.out.validated_file,
        params.format,
        params.hmm_states,
        params.seed,
        params.compute_ld,
        params.max_ld_pairs,
        params.phase_all,
        params.min_maf
    )
    
    // Generate comprehensive report
    generateAnalysisReport(
        runGenomicAlgorithms.out.results,
        runGenomicAlgorithms.out.log,
        runGenomicAlgorithms.out.metrics,
        validateGenomicInput.out.report,
        validateGenomicInput.out.summary
    )
    
    // Final output channels
    runGenomicAlgorithms.out.results.view { "Results: $it" }
    generateAnalysisReport.out.html_report.view { "Report: $it" }
    generateAnalysisReport.out.summary.view { "Summary: $it" }
}
