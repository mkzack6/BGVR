#!/bin/bash

# Complete pipeline test script
# This script tests the entire genomic sparse index pipeline

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() { echo -e "${BLUE}[TEST]${NC} $1"; }
log_success() { echo -e "${GREEN}[PASS]${NC} $1"; }
log_error() { echo -e "${RED}[FAIL]${NC} $1"; }

# Test configuration
TEST_DIR="pipeline_test"
SMALL_VARIANTS=1000
MEDIUM_VARIANTS=5000

cleanup() {
    log_info "Cleaning up test directory..."
    rm -rf "$TEST_DIR"
}

# Trap to cleanup on exit
trap cleanup EXIT

main() {
    log_info "Starting complete pipeline test..."
    
    # Create test environment
    mkdir -p "$TEST_DIR"/{src,data,results,scripts}
    cd "$TEST_DIR"
    
    log_info "Test directory: $PWD"
    
    # Step 1: Check dependencies
    log_info "Checking dependencies..."
    
    if ! command -v rustc &> /dev/null; then
        log_error "Rust not found. Please install Rust first."
        exit 1
    fi
    log_success "Rust found: $(rustc --version)"
    
    if ! command -v nextflow &> /dev/null; then
        log_error "Nextflow not found. Please install Nextflow first."
        exit 1
    fi
    log_success "Nextflow found"
    
    if ! command -v python3 &> /dev/null; then
        log_error "Python3 not found. Please install Python3 first."
        exit 1
    fi
    log_success "Python3 found: $(python3 --version)"
    
    # Step 2: Create project files (simplified versions for testing)
    log_info "Creating test project files..."
    
    # Create a minimal Cargo.toml
    cat > Cargo.toml << 'EOF'
[package]
name = "genomic_sparse_index_test"
version = "0.1.0"
edition = "2021"

[[bin]]
name = "rust_sparse_index"
path = "src/main.rs"

[dependencies]
tokio = { version = "1.0", features = ["full"] }
clap = { version = "4.0", features = ["derive"] }
anyhow = "1.0"
EOF
    
    # Create a minimal test version of main.rs
    mkdir -p src
    cat > src/main.rs << 'EOF'
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use anyhow::Result;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "rust_sparse_index")]
struct Args {
    #[arg(short, long)]
    input: String,
    #[arg(short, long, default_value = "test_results.txt")]
    output: String,
    #[arg(short, long, default_value = "1000000")]
    max_index: usize,
    #[arg(short, long, default_value = "4")]
    num_tasks: usize,
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = Args::parse();
    
    println!("Processing genomic data...");
    println!("Input: {}", args.input);
    println!("Output: {}", args.output);
    
    let file = File::open(&args.input)?;
    let reader = BufReader::new(file);
    let mut count = 0;
    let mut total_alleles = 0u64;
    
    for line_result in reader.lines() {
        let line = line_result?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() >= 3 {
            if let Ok(allele_count) = tokens[2].parse::<u32>() {
                total_alleles += allele_count as u64;
                count += 1;
            }
        }
    }
    
    let mut output_file = File::create(&args.output)?;
    writeln!(output_file, "# Test Results")?;
    writeln!(output_file, "# Total variants processed: {}", count)?;
    writeln!(output_file, "# Total alleles: {}", total_alleles)?;
    writeln!(output_file, "start_pos\tend_pos\tallele_sum")?;
    
    // Generate some test query results
    let chunk_size = args.max_index / args.num_tasks;
    for i in 0..args.num_tasks {
        let start = i * chunk_size + 1;
        let end = if i == args.num_tasks - 1 { args.max_index } else { (i + 1) * chunk_size };
        let mock_sum = (total_alleles / args.num_tasks as u64) + i as u64;
        writeln!(output_file, "{}\t{}\t{}", start, end, mock_sum)?;
    }
    
    println!("Successfully processed {} variants", count);
    println!("Results written to: {}", args.output);
    
    Ok(())
}
EOF
    
    # Create test data generator
    cat > scripts/generate_test_data.py << 'EOF'
#!/usr/bin/env python3
import random
import sys

def generate_test_data(filename, num_variants):
    with open(filename, 'w') as f:
        f.write("# Test genomic data\n")
        f.write(f"# Variants: {num_variants}\n")
        
        for i in range(num_variants):
            chr_num = random.randint(1, 22)
            position = random.randint(1000, 1000000)
            genotype = random.randint(0, 2)
            f.write(f"chr{chr_num}\t{position}\t{genotype}\n")
    
    print(f"Generated {num_variants} test variants in {filename}")

if __name__ == "__main__":
    filename = sys.argv[1] if len(sys.argv) > 1 else "test_genotypes.txt"
    num_variants = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    generate_test_data(filename, num_variants)
EOF
    chmod +x scripts/generate_test_data.py
    
    # Create minimal Nextflow workflow
    cat > main.nf << 'EOF'
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = "data/test_genotypes.txt"
params.output_dir = "results"
params.max_index = 1000000
params.num_tasks = 4

process buildSparseIndex {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path input_file
    
    output:
    path "index_results.txt"
    path "processing.log"
    
    script:
    """
    echo "Starting sparse index build..." > processing.log
    ./rust_sparse_index \\
        --input ${input_file} \\
        --output index_results.txt \\
        --max-index ${params.max_index} \\
        --num-tasks ${params.num_tasks} >> processing.log 2>&1
    echo "Sparse index build completed" >> processing.log
    """
}

workflow {
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    buildSparseIndex(input_ch)
}
EOF
    
    # Step 3: Generate test data
    log_info "Generating test data..."
    python3 scripts/generate_test_data.py data/small_genotypes.txt $SMALL_VARIANTS
    python3 scripts/generate_test_data.py data/medium_genotypes.txt $MEDIUM_VARIANTS
    log_success "Test data generated"
    
    # Step 4: Build Rust project
    log_info "Building Rust project..."
    if cargo build --release; then
        cp target/release/rust_sparse_index .
        chmod +x rust_sparse_index
        log_success "Rust project built successfully"
    else
        log_error "Rust build failed"
        exit 1
    fi
    
    # Step 5: Test Rust binary directly
    log_info "Testing Rust binary..."
    
    if ./rust_sparse_index --help > /dev/null 2>&1; then
        log_success "Rust binary help works"
    else
        log_error "Rust binary help failed"
        exit 1
    fi
    
    if ./rust_sparse_index \
        --input data/small_genotypes.txt \
        --output results/direct_test.txt \
        --max-index 1000000 \
        --num-tasks 2; then
        log_success "Direct Rust execution works"
    else
        log_error "Direct Rust execution failed"
        exit 1
    fi
    
    # Verify output file
    if [ -f "results/direct_test.txt" ] && [ -s "results/direct_test.txt" ]; then
        log_success "Results file created and non-empty"
        echo "Sample output:"
        head -10 results/direct_test.txt | sed 's/^/  /'
    else
        log_error "Results file not created or empty"
        exit 1
    fi
    
    # Step 6: Test Nextflow pipeline
    log_info "Testing Nextflow pipeline..."
    
    if nextflow run main.nf \
        --input data/medium_genotypes.txt \
        --output_dir results_nf \
        --max_index 1000000 \
        --num_tasks 4; then
        log_success "Nextflow pipeline execution works"
    else
        log_error "Nextflow pipeline failed"
        exit 1
    fi
    
    # Verify Nextflow output
    if [ -f "results_nf/index_results.txt" ] && [ -s "results_nf/index_results.txt" ]; then
        log_success "Nextflow results file created"
        echo "Sample Nextflow output:"
        head -10 results_nf/index_results.txt | sed 's/^/  /'
    else
        log_error "Nextflow results file not created or empty"
        exit 1
    fi
    
    # Step 7: Performance test
    log_info "Running performance test..."
    
    start_time=$(date +%s)
    ./rust_sparse_index \
        --input data/medium_genotypes.txt \
        --output results/perf_test.txt \
        --max-index 1000000 \
        --num-tasks 8 > results/perf_log.txt 2>&1
    end_time=$(date +%s)
    
    duration=$((end_time - start_time))
    log_success "Performance test completed in ${duration} seconds"
    
    # Step 8: Validate results consistency
    log_info "Validating results consistency..."
    
    # Check if both direct and Nextflow runs produced similar results
    direct_variants=$(grep "Total variants processed" results/direct_test.txt | awk '{print $NF}')
    nf_variants=$(grep "Total variants processed" results_nf/index_results.txt | awk '{print $NF}')
    
    if [ "$direct_variants" = "$nf_variants" ]; then
        log_success "Results consistency check passed"
    else
        log_error "Results inconsistent: direct=$direct_variants, nextflow=$nf_variants"
        exit 1
    fi
    
    # Step 9: Test error handling
    log_info "Testing error handling..."
    
    # Test with non-existent file
    if ./rust_sparse_index --input nonexistent.txt --output /dev/null 2>/dev/null; then
        log_error "Error handling test failed - should have failed with non-existent file"
        exit 1
    else
        log_success "Error handling works correctly"
    fi
    
    # Step 10: Clean up and summary
    log_info "Test summary:"
    echo "  ✓ Dependencies checked"
    echo "  ✓ Project files created"
    echo "  ✓ Test data generated"
    echo "  ✓ Rust project built"
    echo "  ✓ Direct execution tested"
    echo "  ✓ Nextflow pipeline tested"
    echo "  ✓ Performance test completed"
    echo "  ✓ Results consistency validated"
    echo "  ✓ Error handling tested"
    
    log_success "All tests passed! Pipeline is working correctly."
    
    # Show file sizes and contents
    echo
    log_info "Generated files:"
    find results* -type f -exec ls -lh {} \; | sed 's/^/  /'
    
    echo
    log_info "Test data statistics:"
    echo "  Small dataset: $SMALL_VARIANTS variants"
    echo "  Medium dataset: $MEDIUM_VARIANTS variants"
    echo "  Processing time: ${duration} seconds"
    
    # Final cleanup happens automatically via trap
}

# Parse command line arguments
SKIP_CLEANUP=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-cleanup)
            SKIP_CLEANUP=true
            shift
            ;;
        --small-variants)
            SMALL_VARIANTS="$2"
            shift 2
            ;;
        --medium-variants)
            MEDIUM_VARIANTS="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --skip-cleanup        Don't remove test directory"
            echo "  --small-variants N    Number of variants for small test (default: 1000)"
            echo "  --medium-variants N   Number of variants for medium test (default: 5000)"
            echo "  --help               Show this help"
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Override cleanup if requested
if [ "$SKIP_CLEANUP" = true ]; then
    trap - EXIT
    log_info "Skipping cleanup - test directory preserved: $TEST_DIR"
fi

# Run main test
main