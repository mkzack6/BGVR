### **2.5. Introduction to Machine Learning in Bioinformatics**  
**experiment_2_5_1**  

This Rust program demonstrates the application of a deep neural network using the `tch-rs` crate (Rust bindings for PyTorch) to classify genomic or transcriptomic data. The example assumes a CSV file containing gene expression values, with sample labels indicating the presence or absence of a particular disease phenotype. The goal is to train a multi-layer perceptron (MLP) to distinguish between positive and negative cases based on expression profiles.  

The dataset, **"gene_expression.csv"**, consists of numerical values (one row per sample) with a final column as the classification label (0 or 1). The program reads the CSV, constructs an MLP, and trains a binary classifier using mini-batch processing. While this example focuses on a single dataset, real-world applications often involve large-scale multi-omics data requiring parallel or distributed computing for scalability. Nonetheless, this implementation highlights the fundamental steps in using `tch-rs` for bioinformatics.  

### **Files and Contents**  
- **main.rs** (Rust script)  
- **Python Code for Synthesizing gene_expression.csv.ipynb** (Python script)  
- **gene_expression.csv** ([Google Drive link])  
- **Cargo.toml** (Rust dependencies configuration)  
- **output.txt** (Training output log)  

### **How to Run**  
Execute the following command to run the program and save the output:  
```sh
cargo run | tee output.txt
```  

### **Dependencies**  
```toml
[dependencies]  
tch = "0.19.0"  
```  

---

## **Explanation of the Output**  

### **1. Training Process**  
The neural network is a fully connected model (FCNN) with:  
- **Input Layer**: 10,000 features (gene expression values)  
- **Hidden Layers**: 128 and 64 neurons with ReLU activation  
- **Output Layer**: 1 neuron for binary classification  
- **Loss Function**: Binary Cross-Entropy with Logits (BCEWithLogitsLoss)  
- **Optimizer**: Adam (learning rate = 1e-3)  
- **Batch Size**: 64 samples  

### **2. Validation Loss Progression**  
| Epoch | Validation Loss |  
|-------|----------------|  
| 1     | 71.8563        |  
| 2     | 59.6753        |  
| 3     | 53.6821        |  
| 4     | 51.4969        |  
| 5     | 46.6416        |  
| 6     | 47.4496        |  
| 7     | 38.8133        |  
| 8     | 36.7464        |  
| 9     | 31.4197        |  
| 10    | 25.9429        |  

The loss decreases from 71.8563 in Epoch 1 to 25.9429 in Epoch 10, indicating that the model is learning effectively. The slight increase at Epoch 6 (47.4496) may be due to noisy data, an unstable learning rate, or a challenging batch of samples. However, the downward trend resumes afterward.  

### **3. Key Insights**  
- **Effective Learning**: The steady decline in validation loss suggests successful model training.  
- **Further Optimization**: Additional training epochs might further improve performance.  
- **Overfitting Risk**: If the loss stagnates or increases over time, techniques like dropout or early stopping may be necessary.  
- **Final Evaluation**: The model should be tested on an independent dataset to ensure generalizability.  

### **Conclusion**  
The trained neural network successfully learns from the gene expression dataset, significantly reducing validation loss over 10 epochs. This indicates the model is effectively capturing meaningful patterns in the data and improving its predictions.
