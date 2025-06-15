#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Population Genomics Analysis Pipeline
 * ====================================
 * 
 * This pipeline performs comprehensive population genomics analysis including:
 * - Data generation (synthetic)
 * - Quality control filtering
 * - Principal Component Analysis (PCA)
 * - Population structure analysis
 * - Hardy-Weinberg equilibrium testing
 * - Fst calculations
 */

// Default parameters
params.output_dir = "results"
params.data_dir = "data"
params.individuals = 1000
params.snps = 10000
params.populations = 3
params.maf_threshold = 0.05
params.missing_threshold = 0.1
params.sample_missing_threshold = 0.1
params.num_pcs = 10
params.hwe_test = true
params.population_structure = true
params.help = false

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

def helpMessage() {
    log.info"""
    Population Genomics Analysis Pipeline
    ====================================
    
    Usage:
    nextflow run main.nf [options]
    
    Options:
    --output_dir                Output directory (default: ${params.output_dir})
    --data_dir                  Data directory (default: ${params.data_dir})
    --individuals               Number of individuals to generate (default: ${params.individuals})
    --snps                      Number of SNPs to generate (default: ${params.snps})
    --populations               Number of populations (default: ${params.populations})
    --maf_threshold             Minor allele frequency threshold (default: ${params.maf_threshold})
    --missing_threshold         Missing data threshold per variant (default: ${params.missing_threshold})
    --sample_missing_threshold  Missing data threshold per sample (default: ${params.sample_missing_threshold})
    --num_pcs                   Number of principal components (default: ${params.num_pcs})
    --hwe_test                  Perform Hardy-Weinberg tests (default: ${params.hwe_test})
    --population_structure      Perform population structure analysis (default: ${params.population_structure})
    --help                      Show this help message
    
    Example:
    nextflow run main.nf --individuals 2000 --snps 20000 --populations 5
    """.stripIndent()
}

// Log pipeline parameters
log.info """
Population Genomics Analysis Pipeline
====================================
Output directory      : ${params.output_dir}
Data directory         : ${params.data_dir}
Number of individuals  : ${params.individuals}
Number of SNPs         : ${params.snps}
Number of populations  : ${params.populations}
MAF threshold          : ${params.maf_threshold}
Missing threshold      : ${params.missing_threshold}
Sample missing thresh. : ${params.sample_missing_threshold}
Number of PCs          : ${params.num_pcs}
HWE testing           : ${params.hwe_test}
Population structure  : ${params.population_structure}
"""

/*
 * Process to generate synthetic genomic data
 */
process GENERATE_DATA {
    publishDir params.data_dir, mode: 'copy'
    
    input:
    val individuals
    val snps
    val populations
    
    output:
    path "genotypes.csv", emit: genotypes
    path "sample_metadata.csv", emit: metadata
    path "variant_metadata.csv", emit: variants
    
    script:
    """
    echo "Generating synthetic genomic data..."
    
    # Check if Rust binary exists, if not compile it
    if [ ! -f "../target/release/generate_data" ]; then
        echo "Compiling Rust data generator..."
        cd ..
        cargo build --release --bin generate_data
        cd -
    fi
    
    # Generate the data
    ../target/release/generate_data \\
        --individuals ${individuals} \\
        --snps ${snps} \\
        --populations ${populations} \\
        --output .
        
    echo "Data generation completed successfully!"
    """
}

/*
 * Process to perform quality control and filtering
 */
process QUALITY_CONTROL {
    publishDir "${params.output_dir}/qc", mode: 'copy'
    
    input:
    path genotypes
    path metadata
    val maf_threshold
    val missing_threshold
    val sample_missing_threshold
    
    output:
    path "qc_report.txt", emit: qc_report
    path "filtered_genotypes.csv", emit: filtered_genotypes
    path "qc_summary.log", emit: qc_log
    
    script:
    """
    echo "Performing quality control filtering..." > qc_summary.log
    
    # Check if Rust binary exists, if not compile it
    if [ ! -f "../target/release/population_analysis" ]; then
        echo "Compiling Rust analysis binary..." >> qc_summary.log
        cd ..
        cargo build --release --bin population_analysis
        cd -
    fi
    
    # Run QC-only analysis
    ../target/release/population_analysis \\
        --genotypes ${genotypes} \\
        --metadata ${metadata} \\
        --output . \\
        --maf-threshold ${maf_threshold} \\
        --missing-threshold ${missing_threshold} \\
        --sample-missing-threshold ${sample_missing_threshold} \\
        --num-pcs 5 2>&1 | tee -a qc_summary.log
        
    # Copy filtered data for next step
    cp ${genotypes} filtered_genotypes.csv
    
    echo "Quality control completed!" >> qc_summary.log
    """
}

/*
 * Process to perform Principal Component Analysis
 */
process PCA_ANALYSIS {
    publishDir "${params.output_dir}/pca", mode: 'copy'
    
    input:
    path genotypes
    path metadata
    val num_pcs
    val maf_threshold
    val missing_threshold
    val sample_missing_threshold
    
    output:
    path "pca_eigenvalues.csv", emit: eigenvalues
    path "pca_scores.csv", emit: scores
    path "pca_analysis.log", emit: pca_log
    
    script:
    """
    echo "Performing Principal Component Analysis..." > pca_analysis.log
    
    ../target/release/population_analysis \\
        --genotypes ${genotypes} \\
        --metadata ${metadata} \\
        --output . \\
        --maf-threshold ${maf_threshold} \\
        --missing-threshold ${missing_threshold} \\
        --sample-missing-threshold ${sample_missing_threshold} \\
        --num-pcs ${num_pcs} 2>&1 | tee -a pca_analysis.log
        
    echo "PCA analysis completed!" >> pca_analysis.log
    """
}

/*
 * Process to perform population structure analysis
 */
process POPULATION_STRUCTURE {
    publishDir "${params.output_dir}/population_structure", mode: 'copy'
    
    input:
    path genotypes
    path metadata
    val maf_threshold
    val missing_threshold
    val sample_missing_threshold
    
    output:
    path "population_allele_frequencies.csv", emit: allele_freqs
    path "fst_matrix.csv", emit: fst_matrix
    path "population_structure.log", emit: pop_log
    
    when:
    params.population_structure
    
    script:
    """
    echo "Performing population structure analysis..." > population_structure.log
    
    ../target/release/population_analysis \\
        --genotypes ${genotypes} \\
        --metadata ${metadata} \\
        --output . \\
        --maf-threshold ${maf_threshold} \\
        --missing-threshold ${missing_threshold} \\
        --sample-missing-threshold ${sample_missing_threshold} \\
        --num-pcs 5 \\
        --population-structure 2>&1 | tee -a population_structure.log
        
    echo "Population structure analysis completed!" >> population_structure.log
    """
}

/*
 * Process to perform Hardy-Weinberg equilibrium testing
 */
process HWE_TESTING {
    publishDir "${params.output_dir}/hwe", mode: 'copy'
    
    input:
    path genotypes
    path metadata
    val maf_threshold
    val missing_threshold
    val sample_missing_threshold
    
    output:
    path "hwe_tests.csv", emit: hwe_results
    path "hwe_analysis.log", emit: hwe_log
    
    when:
    params.hwe_test
    
    script:
    """
    echo "Performing Hardy-Weinberg equilibrium testing..." > hwe_analysis.log
    
    ../target/release/population_analysis \\
        --genotypes ${genotypes} \\
        --metadata ${metadata} \\
        --output . \\
        --maf-threshold ${maf_threshold} \\
        --missing-threshold ${missing_threshold} \\
        --sample-missing-threshold ${sample_missing_threshold} \\
        --num-pcs 5 \\
        --hwe-test 2>&1 | tee -a hwe_analysis.log
        
    echo "HWE testing completed!" >> hwe_analysis.log
    """
}

/*
 * Process to generate summary report
 */
process GENERATE_REPORT {
    publishDir params.output_dir, mode: 'copy'
    
    input:
    path qc_report
    path pca_eigenvalues
    path pca_scores
    path allele_freqs
    path fst_matrix
    path hwe_results
    
    output:
    path "final_report.html", emit: report
    path "summary_statistics.txt", emit: summary
    
    script:
    """
    echo "Generating comprehensive analysis report..." > summary_statistics.txt
    echo "=========================================" >> summary_statistics.txt
    echo "" >> summary_statistics.txt
    
    # Add QC summary
    echo "Quality Control Summary:" >> summary_statistics.txt
    if [ -f "${qc_report}" ]; then
        cat ${qc_report} >> summary_statistics.txt
    fi
    echo "" >> summary_statistics.txt
    
    # Add PCA summary
    echo "Principal Component Analysis Summary:" >> summary_statistics.txt
    if [ -f "${pca_eigenvalues}" ]; then
        echo "First 5 principal components:" >> summary_statistics.txt
        head -6 ${pca_eigenvalues} | tail -5 >> summary_statistics.txt
    fi
    echo "" >> summary_statistics.txt
    
    # Add population structure summary
    if [ -f "${fst_matrix}" ]; then
        echo "Population Structure (Fst Matrix):" >> summary_statistics.txt
        cat ${fst_matrix} >> summary_statistics.txt
        echo "" >> summary_statistics.txt
    fi
    
    # Add HWE summary
    if [ -f "${hwe_results}" ]; then
        echo "Hardy-Weinberg Equilibrium Summary:" >> summary_statistics.txt
        echo "Variants significantly deviating from HWE (p < 0.05):" >> summary_statistics.txt
        awk -F',' 'NR>1 && \$4<0.05 {count++} END {print (count ? count : 0) " variants"}' ${hwe_results} >> summary_statistics.txt
        echo "" >> summary_statistics.txt
    fi
    
    # Generate simple HTML report
    cat > final_report.html << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>Population Genomics Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .header { background-color: #f0f0f0; padding: 10px; border-radius: 5px; }
        .section { margin: 20px 0; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <div class="header">
        <h1>Population Genomics Analysis Report</h1>
        <p>Generated on: \$(date)</p>
    </div>
    
    <div class="section">
        <h2>Analysis Overview</h2>
        <p>This report summarizes the results of a comprehensive population genomics analysis pipeline.</p>
    </div>
    
    <div class="section">
        <h2>Files Generated</h2>
        <ul>
            <li>Quality Control Report: qc/qc_report.txt</li>
            <li>PCA Results: pca/pca_eigenvalues.csv, pca/pca_scores.csv</li>
            <li>Population Structure: population_structure/fst_matrix.csv</li>
            <li>Hardy-Weinberg Tests: hwe/hwe_tests.csv</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>Summary Statistics</h2>
        <pre>
EOF
    cat summary_statistics.txt >> final_report.html
    cat >> final_report.html << 'EOF'
        </pre>
    </div>
</body>
</html>
EOF
    
    echo "Report generation completed!"
    """
}

/*
 * Main workflow
 */
workflow {
    // Generate synthetic data
    GENERATE_DATA(
        params.individuals,
        params.snps,
        params.populations
    )
    
    // Quality control
    QUALITY_CONTROL(
        GENERATE_DATA.out.genotypes,
        GENERATE_DATA.out.metadata,
        params.maf_threshold,
        params.missing_threshold,
        params.sample_missing_threshold
    )
    
    // Principal Component Analysis
    PCA_ANALYSIS(
        GENERATE_DATA.out.genotypes,
        GENERATE_DATA.out.metadata,
        params.num_pcs,
        params.maf_threshold,
        params.missing_threshold,
        params.sample_missing_threshold
    )
    
    // Population structure analysis (conditional)
    if (params.population_structure) {
        POPULATION_STRUCTURE(
            GENERATE_DATA.out.genotypes,
            GENERATE_DATA.out.metadata,
            params.maf_threshold,
            params.missing_threshold,
            params.sample_missing_threshold
        )
    }
    
    // Hardy-Weinberg equilibrium testing (conditional)
    if (params.hwe_test) {
        HWE_TESTING(
            GENERATE_DATA.out.genotypes,
            GENERATE_DATA.out.metadata,
            params.maf_threshold,
            params.missing_threshold,
            params.sample_missing_threshold
        )
    }
    
    // Generate final report
    GENERATE_REPORT(
        QUALITY_CONTROL.out.qc_report,
        PCA_ANALYSIS.out.eigenvalues,
        PCA_ANALYSIS.out.scores,
        params.population_structure ? POPULATION_STRUCTURE.out.allele_freqs : file('NO_FILE'),
        params.population_structure ? POPULATION_STRUCTURE.out.fst_matrix : file('NO_FILE'),
        params.hwe_test ? HWE_TESTING.out.hwe_results : file('NO_FILE')
    )
}

/*
 * Workflow completion message
 */
workflow.onComplete {
    println ""
    println "Pipeline completed!"
    println "Results are available in: ${params.output_dir}"
    println ""
    if (workflow.success) {
        println "Analysis completed successfully! 🎉"
        println ""
        println "Key outputs:"
        println "- Quality control report: ${params.output_dir}/qc/qc_report.txt"
        println "- PCA results: ${params.output_dir}/pca/"
        if (params.population_structure) {
            println "- Population structure: ${params.output_dir}/population_structure/"
        }
        if (params.hwe_test) {
            println "- Hardy-Weinberg tests: ${params.output_dir}/hwe/"
        }
        println "- Final report: ${params.output_dir}/final_report.html"
    } else {
        println "Pipeline failed. Check the error messages above."
    }
}
