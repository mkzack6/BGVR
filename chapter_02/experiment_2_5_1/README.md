# 2.5. Genomic Classification
## Experiment 2.5.1: Binary Classification of Gene Expression Data
This project provides a Rust program (`rust_genomics`) to perform binary classification on gene expression data, such as those derived from RNA sequencing or microarray experiments. The program is designed for educational purposes, demonstrating machine learning in bioinformatics using Rust’s high-performance capabilities and the tch crate for tensor operations. It processes gene expression data in CSV format, trains a neural network to predict binary labels (e.g., disease vs. healthy), and outputs training metrics. The implementation leverages Rust’s integration with LibTorch for efficient model training, making it suitable for analyzing genomic datasets in high-performance computing (HPC) environments.
## Breakdown of the Project
### Dependencies
The project consists of a single Rust program (rust_genomics). Dependencies are assumed to be defined in Cargo.toml:
```toml
[dependencies]
tch = "0.19.0"
csv = "1.3"
anyhow = "1.0"
```

- **tch**: Provides tensor operations and neural network training via LibTorch, enabling machine learning on gene expression data.
- **csv**: Supports parsing of CSV files containing gene expression data and labels.
- **anyhow**: Facilitates flexible error handling for robust program execution.

### Feature Flags
The program has no feature flags, using a single workflow with fixed parameters:

- **Default**: Expects a CSV file (`gene_expression.csv`) with 10,000 gene expression features and a binary label per row.
- **Data**: Produces training metrics (logits mean, loss, validation loss) and saves results to `output.txt`.

The program has one main process:

- **rust_genomics**: Loads gene expression data, trains a neural network over 10 epochs, logs training and validation metrics, and saves results.

### Implementation Modules
Rust Program (`src/main.rs`):

- **Data Structures**:
-   **GeneExpression**: Represents a row of gene expression data with 10,000 floating-point features (e.g., `[0.5243, 0.7812, ..., 0.6723]`) and a binary label (0 or 1).
-   **Tensor**: Uses `tch` tensors to store feature matrices and labels for model training (e.g., shape `[num_samples, 10000]` for features).


- **Model Training**:
-   Implements a neural network (assumed to be a simple feedforward network) to predict binary labels.
-   Uses a loss function (likely binary cross-entropy) to optimize model weights.
-   Processes data in batches (e.g., 65 batches per epoch), logging logits mean and loss for batches 0 and 64.
-   Computes validation loss after each epoch to assess generalization.


- **Input/Output**:
-   Reads `gene_expression.csv`, expecting a header (`gene1,gene2,...,gene10000,label`) and rows with 10,001 columns (10,000 features + 1 label).
-   Outputs training metrics to the console (logits mean, loss, validation loss).
-   Saves detailed results to `output.txt` (contents depend on implementation, e.g., final model weights or predictions).


- **Main Function**:
-   Parses the CSV file, converts data to tensors, trains the model for 10 epochs, and saves results.
-   Outputs metrics for ~100 samples (based on observed batch sizes), producing an `output.txt` file of unknown size.



## Program Execution
The program processes gene expression data from `gene_expression.csv` to train a neural network and produce training metrics. Steps to execute:
### Prerequisites

**Install Rust/Cargo**:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

**Install LibTorch**:

Download and extract LibTorch (CPU version) for your system from https://pytorch.org/.
Set the LIBTORCH environment variable:
```bash
export LIBTORCH=/home/zack/libtorch
export LD_LIBRARY_PATH=$LIBTORCH/lib:$LD_LIBRARY_PATH
```

Example for Ubuntu (WSL2):wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.0%2Bcpu.zip
unzip libtorch-cxx11-abi-shared-with-deps-2.0.0+cpu.zip -d /home/zack/


**Install Python** (for CSV generation):

Required to generate gene_expression.csv if not using a real dataset:
```bash
sudo apt update
sudo apt install python3 python3-pip`
```


### Setup

Ensure the `rust_genomics/` directory (containing `Cargo.toml` and `src/main.rs`) is in your project directory:`~/rust_genomics/`

```csv
Verify or create gene_expression.csv in the same directory. Example structure:gene1,gene2,gene3,...,gene9998,gene9999,gene10000,label
0.5243,0.7812,0.3467,...,0.9123,0.4512,0.6723,1
0.2314,0.8976,0.1245,...,0.3345,0.7890,0.5634,0
...
```

Generate `gene_expression.csv` using the following Python script (`generate_gene_expression.py`):
```python
import csv
import random

num_genes = 10000
num_rows = 100

header = [f"gene{i}" for i in range(1, num_genes + 1)] + ["label"]
data = []
for _ in range(num_rows):
    row = [round(random.uniform(0.0, 1.0), 4) for _ in range(num_genes)]
    row.append(random.randint(0, 1))
    data.append(row)

with open("gene_expression.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(data)
```
Run the script in WSL2:
```bash
python3 generate_gene_expression.py
```



### Run Program

Navigate to the project directory:
```bash
cd ~/rust_genomics/
```


Build the program:
```bash
cargo build --release
```

Run the program:
```bash
cargo run --release
```

Note: The program assumes `gene_expression.csv` is in the current directory. If a different path is required, modify the Rust code or pass the path as an argument (if supported).

### Example Output
`First 10 feature values: Tensor[[1, 10000], Float]
Epoch 1, Batch 0: Logits Mean [-0.32026344537734985], Loss [0.6811246871948242]
Epoch 1, Batch 64: Logits Mean [-0.6296641826629639], Loss [0.7382960319519043]
Epoch 1, Validation loss: 1.4194
...
Epoch 10, Batch 64: Logits Mean [-0.1613333523273468], Loss [0.27567702531814575]
Epoch 10, Validation loss: 0.7602
Training complete. Check output file at: /home/zack/rust_genomics/output.txt`

### Generated Files

- **output.txt**: Contains detailed training results (e.g., final model weights, predictions, or metrics). Location: `/home/zack/rust_genomics/output.txt`.
- **gene_expression.csv** (if generated): Input file with 100 rows, 10,001 columns (~10 MB).

## Potential Errors and Solutions

### Incorrect Feature Count in CSV:

- **Error**: `Warning: Skipping row with incorrect feature count: [0.0, 1.0]`
- **Cause**: The CSV file has rows with fewer columns (e.g., 2) than expected (10,001).
- **Solution**: Ensure `gene_expression.csv` has 10,000 gene features + 1 label per row. Regenerate the CSV using the provided Python script or verify the file’s structure:
```bash
head -n 5 gene_expression.csv
```

- **Torch Stack Error**:

- **Error**: `thread 'main' panicked at ...: called` Result::unwrap() `on an` Err `value: Torch("stack expects a non-empty TensorList..."`)
- **Cause**: No valid rows were processed (e.g., all rows skipped due to incorrect feature count), resulting in an empty tensor list for the stack operation.
- **Solution**: Fix the CSV file to ensure all rows have the correct number of columns. Verify data parsing in the Rust code (e.g., check CSV parsing logic in `src/main.rs`).


- **LibTorch Not Found**:

- **Error**: `error: could not find native library` libtorch ...
- **Cause**: LibTorch is not installed or the `LIBTORCH` environment variable is not set.
- **Solution**: Install LibTorch and set environment variables:
```bash
export LIBTORCH=/home/zack/libtorch
export LD_LIBRARY_PATH=$LIBTORCH/lib:$LD_LIBRARY_PATH
```

- **Python Script Execution Issues**:

- **Error**: `python3 command not found` or script fails to generate CSV.
- **Cause**: Python is not installed, or the script has syntax/file path errors.
- **Solution**: Install Python and verify the script:
```bash
sudo apt install python3 python3-pip
python3 generate_gene_expression.py
```
Ensure the script is in the correct directory (`~/rust_genomics/`).


- **High Validation Loss or Training Instability**:

- **Error**: Validation loss spikes (e.g., `5.0425` in Epoch 2) or remains high (`0.7602` in Epoch 10).
- **Cause**: Randomly generated data lacks meaningful patterns, small dataset size (100 rows), or suboptimal hyperparameters.
- **Solution**: Use a real gene expression dataset with more rows.
-   **Tune hyperparameters** (e.g., learning rate, batch size, epochs).
-   **Increase dataset size** (e.g., 1,000 rows) by modifying the Python script (`num_rows = 1000`).



## Why is this Project Important?

- **Bioinformatics Workflow**: Demonstrates binary classification of gene expression data, critical for applications like disease diagnosis or biomarker discovery.
- **Machine Learning**: Leverages Rust and LibTorch for high-performance model training, showcasing Rust’s capabilities in scientific computing.
- **Data Processing**: Handles large-scale CSV data, a common format in genomics, with robust error handling.
- **Educational Value**: Teaches students Rust, tensor operations, and machine learning in a bioinformatics context, aligning with Chapter 5’s focus on genomic analysis.
- **Scalability**: Suitable for HPC environments, where large genomic datasets require efficient processing.

## Next Steps / Improvements

- **Real Data Integration**: Replace the random CSV with a real gene expression dataset (e.g., from GEO or TCGA) to improve model performance.
- **Model Tuning**: Adjust learning rate, batch size, or model architecture to reduce validation loss and stabilize training.
- **Metrics Reporting**: Add accuracy, precision, recall, and F1-score to the output for better evaluation.
- **Dataset Size**: Increase the number of samples (e.g., 1,000 or 10,000 rows) to improve training stability.
- **Visualization**: Plot training/validation loss using Python/Matplotlib to visualize progress.
- **Pipeline Integration**: Wrap in a Nextflow or Snakemake pipeline for automated HPC workflows.
- **Output Formats**: Support additional output formats (e.g., CSV or JSON) for compatibility with downstream analysis tools.

