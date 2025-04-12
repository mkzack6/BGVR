# 3.1 Introduction to Data Structures and Algorithms

**Experiment 3_1_4: Conditional Hardware Acceleration in Rust**

This Rust project demonstrates **conditional computation offloading** to different hardware accelerators (CPU, GPU, FPGA) using Rust's feature flags. It provides a unified interface while enabling specialized code paths for different hardware targets at compile time.

## Breakdown of the Project

### 1. Dependencies (Cargo.toml)
The project uses:
* `rand`: A random number generation library.
* `rayon`: A parallel computing library for CPU parallelism.

```toml
[dependencies]
rand = "0.9.0"
rayon = "1.10.0"
[features]
gpu = []
fpga = []
```

### 2. Feature Flags
The project uses Rust's feature flags to enable different hardware targets:
* Default (no flags): Uses CPU implementation
* `gpu`: Enables GPU acceleration path
* `fpga`: Enables FPGA acceleration path

### 3. Hardware-Specific Implementation Modules
The implementation:
1. **GPU Module** (`#[cfg(feature = "gpu")]`):
   * Contains GPU-specific initialization and acceleration functions.
   * In a real implementation, would interface with CUDA, OpenCL, or Vulkan.
   * Example shows a simulated GPU computation doubling input values.

2. **FPGA Module** (`#[cfg(feature = "fpga")]`):
   * Contains FPGA-specific initialization and acceleration functions.
   * In a real implementation, would interface with vendor-specific APIs.
   * Example shows a simulated FPGA computation adding 10.0 to input values.

3. **CPU Fallback** (`#[cfg(not(any(feature = "gpu", feature = "fpga")))]`):
   * Provides a CPU-based implementation when no accelerators are enabled.
   * Uses parallel processing via Rayon if needed.
   * Example shows multiplication of input values by 1.1.

Example: Processing a vector `[1.0, 2.0, 3.0, 4.0, 5.0]` with different targets:
```
GPU: [2.0, 4.0, 6.0, 8.0, 10.0]  // Each value multiplied by 2.0
FPGA: [11.0, 12.0, 13.0, 14.0, 15.0]  // Each value increased by 10.0
CPU: [1.1, 2.2, 3.3, 4.4, 5.5]  // Each value multiplied by 1.1
```

### 4. Program Execution
* The `main()` function demonstrates how to use conditional compilation.
* Based on enabled features, it calls the appropriate acceleration function.
* Can be built with different targets using cargo features:
  * `cargo build` (CPU only)
  * `cargo build --features gpu` (GPU acceleration)
  * `cargo build --features fpga` (FPGA acceleration)

### Example Output:
```
No GPU/FPGA features enabled. Using CPU fallback in parallel if desired.
CPU-only fallback result: [1.1, 2.2, 3.3000002, 4.4, 5.5]
```

* Output shows the computation was performed using the CPU fallback path.
* The small precision difference in `3.3000002` is typical of floating-point arithmetic.

## Why is this Project Important?
* **Allows a single codebase** to adapt to different hardware environments.
* **Reduces maintenance complexity** by avoiding runtime branching.
* Uses Rust's **zero-cost abstractions** to eliminate runtime overhead.
* **Compile-time specialization** ensures optimal code paths for each target.
* Creates **smaller binaries** by including only the code needed for the selected features.

## Next Steps / Improvements
* **Real Hardware Integration**: Replace simulated functions with actual GPU/FPGA interfaces.
* **Benchmarking Framework**: Add tools to compare performance across different hardware.
* **Auto-Detection**: Add build scripts to automatically detect available hardware.
* **Runtime Fallback**: Add runtime detection as a safety mechanism.
* **Extended Abstractions**: Create higher-level abstractions for complex algorithms.
