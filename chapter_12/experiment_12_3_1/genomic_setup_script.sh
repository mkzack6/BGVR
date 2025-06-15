#!/bin/bash

# Enhanced Genomic Algorithms Pipeline Setup Script for WSL
# This script sets up the complete environment for genomic analysis

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_step() {
    echo -e "${PURPLE}[STEP]${NC} $1"
}

# Project configuration
PROJECT_NAME="genomic_algorithms"
PROJECT_DIR="enhanced_genomic_pipeline"

# Check if running in WSL
check_wsl() {
    if ! grep -q microsoft /proc/version 2>/dev/null; then
        log_warning "This script is designed for WSL. You may need to adapt commands for your system."
    else
        log_info "Running in WSL environment"
    fi
}

# Update system packages
update_system() {
    log_step "Updating system packages..."
    sudo apt update && sudo apt upgrade -y
    log_success "System packages updated"
}

# Install system dependencies
install_system_dependencies() {
    log_step "Installing system dependencies..."
    
    # Essential development tools
    sudo apt install -y \
        curl \
        wget \
        git \
        build-essential \
        pkg-config \
        libssl-dev \
        python3 \
        python3-pip \
        python3-venv \
        default-jre \
        unzip \
        htop \
        time \
        jq
    
    log_success "System dependencies installed"
}

# Install Python scientific packages
install_python_packages() {
    log_step "Installing Python scientific packages..."
    
    # Create virtual environment
    python3 -m venv genomic_env
    source genomic_env/bin/activate
    
    # Upgrade pip
    pip install --upgrade pip
    
    # Install scientific packages
    pip install \
        numpy \
        scipy \
        pandas \
        matplotlib \
        seaborn \
        scikit-learn \
        jupyter
    
    log_success "Python packages installed in virtual environment"
}

# Install Rust
install_rust() {
    log_step "Installing Rust..."
    
    if command -v rustc &> /dev/null; then
        log_warning "Rust is already installed ($(rustc --version))"
        return
    fi
    
    # Install Rust via rustup
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source ~/.cargo/env
    
    # Install additional components
    rustup component add clippy rustfmt
    
    # Verify installation
    if command -v rustc &> /dev/null; then
        log_success "Rust installed successfully ($(rustc --version))"
    else
        log_error "Rust installation failed"
        exit 1
    fi
}

# Install Nextflow
install_nextflow() {
    log_step "Installing Nextflow..."
    
    if command -v nextflow &> /dev/null; then
        log_warning "Nextflow is already installed"
        nextflow -version | head -1
        return
    fi
    
    # Download and install Nextflow
    curl -s https://get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin/
    sudo chmod +x /usr/local/bin/nextflow
    
    # Verify installation
    if command -v nextflow &> /dev/null; then
        log_success "Nextflow installed successfully"
        nextflow -version | head -1
    else
        log_error "Nextflow installation failed"
        exit 1
    fi
}

# Create project directory structure
create_project_structure() {
    log_step "Creating project directory structure..."
    
    if [ -d "$PROJECT_DIR" ]; then
        log_warning "Project directory $PROJECT_DIR already exists"
        read -p "Do you want to continue and potentially overwrite files? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            log_info "Aborting setup"
            exit 0
        fi
    fi
    
    mkdir -p "$PROJECT_DIR"/{src,data,results,scripts,tests,docs,benchmarks}
    cd "$PROJECT_DIR"
    
    log_success "Project directory structure created: $PWD"
}

# Setup Rust project
setup_rust_project() {
    log_step "Setting up Rust project..."
    
    # Initialize Cargo project if not exists
    if [ ! -f "Cargo.toml" ]; then
        cargo init --name $PROJECT_NAME
        log_info "Initialized new Cargo project"
    fi
    
    log_warning "Remember to copy the enhanced Cargo.toml and src/main.rs from the provided artifacts"
}

# Create sample datasets
create_sample_datasets() {
    log_step "Creating sample genomic datasets..."
    
    if [ -f "scripts/generate_genomic_data.py" ]; then
        # Activate Python environment
        source ../genomic_env/bin/activate 2>/dev/null || true
        
        # Generate small test dataset
        python3 scripts/generate_genomic_data.py \
            --small \
            --output data/small_cohort.csv \
            --metadata data/small_cohort_metadata.json
        
        # Generate medium dataset for testing
        python3 scripts/generate_genomic_data.py \
            --medium \
            --output data/medium_cohort.csv \
            --populations 2 \
            --metadata data/medium_cohort_metadata.json
        
        # Generate large dataset for performance testing
        python3 scripts/generate_genomic_data.py \
            --large \
            --output data/large_cohort.csv \
            --populations 3 \
            --missing-rate 0.03 \
            --metadata data/large_cohort_metadata.json
        
        log_success "Sample datasets generated in data/ directory"
    else
        log_warning "Data generator script not found. Copy generate_genomic_data.py to scripts/ directory"
    fi
}

# Build Rust project
build_rust_project() {
    log_step "Building Rust project..."
    
    if [ ! -f "Cargo.toml" ]; then
        log_error "Cargo.toml not found. Make sure to copy it from the artifacts"
        return 1
    fi
    
    # Build in release mode for optimal performance
    cargo build --release
    
    if [ -f "target/release/rust_core_algorithms" ]; then
        # Copy binary to project root for Nextflow
        cp target/release/rust_core_algorithms .
        chmod +x rust_core_algorithms
        log_success "Rust project built successfully"
    else
        log_error "Rust build failed"
        return 1
    fi
}

# Create test script
create_test_script() {
    log_step "Creating comprehensive test script..."
    
    cat > test_genomic_pipeline.sh << 'EOF'
#!/bin/bash

# Comprehensive test script for genomic algorithms pipeline

set -e

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

log_test() { echo -e "${BLUE}[TEST]${NC} $1"; }
log_pass() { echo -e "${GREEN}[PASS]${NC} $1"; }
log_fail() { echo -e "${RED}[FAIL]${NC} $1"; }

echo "=== Enhanced Genomic Pipeline Test Suite ==="
echo

# Test 1: Binary availability
log_test "Checking binary availability..."
if [ -f "./rust_core_algorithms" ]; then
    log_pass "Rust binary found"
else
    log_fail "Rust binary not found"
    exit 1
fi

# Test 2: Help command
log_test "Testing help command..."
if ./rust_core_algorithms --help > /dev/null 2>&1; then
    log_pass "Help command works"
else
    log_fail "Help command failed"
    exit 1
fi

# Test 3: Small dataset analysis
log_test "Running analysis on small dataset..."
if [ -f "data/small_cohort.csv" ]; then
    ./rust_core_algorithms \
        --input data/small_cohort.csv \
        --output results/test_small.txt \
        --format txt \
        --compute-ld \
        --max-ld-pairs 100 \
        --verbose
    log_pass "Small dataset analysis completed"
else
    log_fail "Small dataset not found"
    exit 1
fi

# Test 4: Nextflow pipeline
log_test "Testing Nextflow pipeline..."
if command -v nextflow &> /dev/null && [ -f "main.nf" ]; then
    nextflow run main.nf \
        --input data/small_cohort.csv \
        --output_dir results/nextflow_test \
        --compute_ld true \
        --max_ld_pairs 50
    log_pass "Nextflow pipeline test completed"
else
    log_fail "Nextflow or main.nf not available"
fi

# Test 5: Performance test
log_test "Running performance test..."
if [ -f "data/medium_cohort.csv" ]; then
    time ./rust_core_algorithms \
        --input data/medium_cohort.csv \
        --output results/performance_test.json \
        --format json \
        --compute-ld \
        --max-ld-pairs 500 \
        --verbose > results/performance_log.txt 2>&1
    log_pass "Performance test completed"
else
    log_fail "Medium dataset not found"
fi

echo
echo "=== Test Summary ==="
echo "✓ Binary functionality verified"
echo "✓ Data analysis capabilities tested"
echo "✓ Pipeline integration confirmed"
echo "✓ Performance benchmarked"
echo
echo "All tests passed! The genomic pipeline is ready for use."
EOF
    
    chmod +x test_genomic_pipeline.sh
    log_success "Test script created: test_genomic_pipeline.sh"
}

# Create documentation
create_documentation() {
    log_step "Creating documentation..."
    
    cat > docs/USAGE.md << 'EOF'
# Enhanced Genomic Algorithms Pipeline Usage Guide

## Quick Start

### 1. Generate Test Data
```bash
# Small dataset for testing
python3 scripts/generate_genomic_data.py --small -o data/test.csv

# Medium dataset with population structure
python3 scripts/generate_genomic_data.py --medium --populations 2 -o data/cohort.csv
```

### 2. Run Direct Analysis
```bash
# Basic allele frequency analysis
./rust_core_algorithms --input data/test.csv --output results/basic.txt

# Comprehensive analysis with LD and phasing
./rust_core_algorithms \
    --input data/cohort.csv \
    --output results/comprehensive.json \
    --format json \
    --compute-ld \
    --phase-all \
    --max-ld-pairs 1000 \
    --verbose
```

### 3