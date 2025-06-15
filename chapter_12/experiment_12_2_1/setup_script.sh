#!/bin/bash

# Genomic Sparse Index Pipeline Setup Script for WSL
# This script sets up the complete environment and dependencies

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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

# Check if running in WSL
check_wsl() {
    if ! grep -q microsoft /proc/version; then
        log_warning "This script is designed for WSL. You may need to adapt commands for your system."
    else
        log_info "Running in WSL environment"
    fi
}

# Update system packages
update_system() {
    log_info "Updating system packages..."
    sudo apt update && sudo apt upgrade -y
    log_success "System packages updated"
}

# Install basic dependencies
install_dependencies() {
    log_info "Installing basic dependencies..."
    
    # Essential tools
    sudo apt install -y \
        curl \
        wget \
        git \
        build-essential \
        pkg-config \
        libssl-dev \
        python3 \
        python3-pip \
        default-jre \
        unzip
    
    log_success "Basic dependencies installed"
}

# Install Rust
install_rust() {
    log_info "Installing Rust..."
    
    if command -v rustc &> /dev/null; then
        log_warning "Rust is already installed ($(rustc --version))"
        return
    fi
    
    # Install Rust via rustup
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source ~/.cargo/env
    
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
    log_info "Installing Nextflow..."
    
    if command -v nextflow &> /dev/null; then
        log_warning "Nextflow is already installed ($(nextflow -version | head -1))"
        return
    fi
    
    # Download and install Nextflow
    curl -s https://get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin/
    sudo chmod +x /usr/local/bin/nextflow
    
    # Verify installation
    if command -v nextflow &> /dev/null; then
        log_success "Nextflow installed successfully"
    else
        log_error "Nextflow installation failed"
        exit 1
    fi
}

# Create project directory structure
create_project_structure() {
    log_info "Creating project directory structure..."
    
    PROJECT_DIR="genomic_sparse_index"
    
    if [ -d "$PROJECT_DIR" ]; then
        log_warning "Project directory $PROJECT_DIR already exists"
        read -p "Do you want to continue and potentially overwrite files? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            log_info "Aborting setup"
            exit 0
        fi
    fi
    
    mkdir -p "$PROJECT_DIR"/{src,data,results,scripts}
    cd "$PROJECT_DIR"
    
    log_success "Project directory structure created: $PWD"
}

# Setup Rust project
setup_rust_project() {
    log_info "Setting up Rust project..."
    
    # Initialize Cargo project if not exists
    if [ ! -f "Cargo.toml" ]; then
        cargo init --name genomic_sparse_index
        log_info "Initialized new Cargo project"
    fi
    
    # Note: The Cargo.toml and src/main.rs files should be created from the artifacts
    log_warning "Remember to copy the Cargo.toml and src/main.rs from the provided artifacts"
}

# Generate sample data
generate_sample_data() {
    log_info "Generating sample genomic data..."
    
    # Note: The Python script should be copied from artifacts
    if [ -f "scripts/generate_sample_data.py" ]; then
        python3 scripts/generate_sample_data.py -o data/genotypes.txt -n 50000 -m 5000000
        log_success "Sample data generated: data/genotypes.txt"
    else
        log_warning "Sample data generator not found. Copy generate_sample_data.py to scripts/ directory"
    fi
}

# Build Rust project
build_rust_project() {
    log_info "Building Rust project..."
    
    if [ ! -f "Cargo.toml" ]; then
        log_error "Cargo.toml not found. Make sure to copy it from the artifacts"
        return 1
    fi
    
    # Build in release mode for better performance
    cargo build --release
    
    if [ -f "target/release/rust_sparse_index" ]; then
        # Copy binary to project root for Nextflow
        cp target/release/rust_sparse_index .
        chmod +x rust_sparse_index
        log_success "Rust project built successfully"
    else
        log_error "Rust build failed"
        return 1
    fi
}

# Test the setup
test_setup() {
    log_info "Testing the setup..."
    
    # Test Rust binary
    if [ -f "rust_sparse_index" ]; then
        ./rust_sparse_index --help > /dev/null 2>&1 && \
        log_success "Rust binary works correctly" || \
        log_error "Rust binary test failed"
    else
        log_error "Rust binary not found"
    fi
    
    # Test Nextflow
    if command -v nextflow &> /dev/null; then
        nextflow -version > /dev/null 2>&1 && \
        log_success "Nextflow works correctly" || \
        log_error "Nextflow test failed"
    else
        log_error "Nextflow not found"
    fi
    
    # Check for required files
    local required_files=("Cargo.toml" "src/main.rs" "main.nf")
    for file in "${required_files[@]}"; do
        if [ -f "$file" ]; then
            log_success "Required file found: $file"
        else
            log_warning "Required file missing: $file (copy from artifacts)"
        fi
    done
}

# Display next steps
show_next_steps() {
    log_info "Setup completed! Next steps:"
    echo
    echo "1. Copy the following files from the artifacts to the project directory:"
    echo "   - Cargo.toml (to project root)"
    echo "   - src/main.rs (to src/ directory)"
    echo "   - main.nf (to project root)"
    echo "   - scripts/generate_sample_data.py (to scripts/ directory)"
    echo
    echo "2. Generate sample data:"
    echo "   python3 scripts/generate_sample_data.py -o data/genotypes.txt -n 50000"
    echo
    echo "3. Build the Rust project:"
    echo "   cargo build --release"
    echo "   cp target/release/rust_sparse_index ."
    echo
    echo "4. Test the Rust binary:"
    echo "   ./rust_sparse_index --help"
    echo "   ./rust_sparse_index -i data/genotypes.txt -o results/test_results.txt"
    echo
    echo "5. Run the Nextflow pipeline:"
    echo "   nextflow run main.nf --input data/genotypes.txt --output_dir results"
    echo
    echo "Project directory: $PWD"
    echo
    log_success "All done! Happy genomic data processing!"
}

# Main execution
main() {
    log_info "Starting Genomic Sparse Index Pipeline setup for WSL..."
    echo
    
    check_wsl
    update_system
    install_dependencies
    install_rust
    install_nextflow
    create_project_structure
    setup_rust_project
    
    # Source cargo environment for current session
    source ~/.cargo/env
    
    test_setup
    show_next_steps
}

# Run main function
main "$@"