## 2.4. Pangenome Graph Theorems
Experiment 2_4_1

This Rust project constructs a **pangenome graph** from multiple synthetic FASTA sequences, representing different haplotypes. A **pangenome graph** captures the diversity among genomes by encoding variations within a single graph structure. Your project extracts **overlapping k-mers** (subsequences of length `k`) from input FASTA files, builds a directed graph where nodes represent k-mers, and edges indicate sequence continuity.

---

## **Breakdown of the Project**
### **1. Dependencies (Cargo.toml)**
The project uses:
- **`bio`**: A Rust bioinformatics library for handling FASTA files.
- **`rayon`**: A parallel computing library to speed up processing.

---

### **2. FASTA Input Files (`src/haplotype1.fasta`, etc.)**
Each FASTA file represents a different haplotype (variation of a genome). Example:
```
>haplotype1_chr1
ACGTACGTACGTACGTACGTACGTACGTTTGGGCCCACGTACGTACGTAAAAC
```
The program reads these files and extracts **k-mers** (length `k=21` in this case).

---

### **3. `build_pangenome_graph()` Function**
This function:
1. **Reads FASTA files in parallel** using Rayon.
2. **Extracts k-mers and builds a directed adjacency list**:
   - Each **k-mer** becomes a **node**.
   - Overlapping **(k-1) suffix-prefix relations** define **edges**.
3. **Merges graphs from all haplotypes** to form a **global pangenome graph**.

Example:
For `k=5`, if a sequence is `ACGTACGT`, the extracted k-mers are:
```
ACGTA → CGTAC → GTACG → TACGT
```
Each arrow represents a directed edge in the graph.

---

### **4. `main.rs` Execution**
- Calls `build_pangenome_graph(k=21, haplotypes)`.
- Constructs the graph and prints:
  - The **total number of nodes** (k-mers).
  - Example **edges** in the graph.

**Example Output (`output.txt`):**
```
Constructed a pangenome graph with 142 nodes.
Node: ACGTTTTTGAAAACCCTGGGA -> ["CGTTTTTGAAAACCCTGGGAC"]
Node: AAACCCTGGGACGTACGTACG -> ["AACCCTGGGACGTACGTACGT"]
Node: ACGTAAAACACGTACGTACGT -> ["CGTAAAACACGTACGTACGTA"]
Node: GTACGTTTTTGAAAACCCTGG -> ["TACGTTTTTGAAAACCCTGGG"]
Node: GTACGTACGTACGTACGTACG -> ["TACGTACGTACGTACGTACGT", "TACGTACGTACGTACGTACGT", ...]
```
- **First column** = A node (a k-mer).
- **Second column** = List of **connected k-mers** (edges in the graph).
- Example:
  ```
  ACGTTTTTGAAAACCCTGGGA -> CGTTTTTGAAAACCCTGGGAC
  ```
  means the sequence **"CGTTTTTGAAAACCCTGGGAC"** follows **"ACGTTTTTGAAAACCCTGGGA"** in the dataset.

---

## **Why is this Project Important?**
- **Pangenome graphs** capture **genetic variation** across multiple genomes.
- Helps in **genome analysis**, **variant detection**, and **comparative genomics**.
- Can be extended for **machine learning on genomic sequences**.

---

## **Next Steps / Improvements**
- **Graph Visualization**: Use a tool like `Graphviz` to visualize nodes and edges.
- **Graph Optimization**: Apply compression techniques for large-scale genomes.
- **Biological Interpretation**: Integrate metadata (e.g., functional annotations).

