# ETO-HEOM (GPU Version)

This directory contains the **GPU-accelerated implementation** of the Extended Thiele Oscillator Hierarchical Equations of Motion (**ETO-HEOM**) framework, optimized for simulating open quantum systems in the context of 2D electronic spectroscopy.

---

## Structure

- **`lib/`**  
  Contains the static CUDA-compiled library (`libgpu2des.a`) and object files.

- **`src/`**  
  Contains the CUDA-based source files for:
  - `GPU_2DES`: executable for 2D electronic spectroscopy simulations
  - `GPU_dynamics`: executable for population dynamics or real-time evolution

---

## GPU & CUDA Requirements

- This version **requires an NVIDIA GPU** and uses:
  - **CUDA Toolkit 11.8**
  - `nvcc` for compilation
  - CUDA libraries for parallel matrix operations and integration kernels

Make sure your system supports CUDA 11.8 and the necessary GPU drivers are correctly installed.

---

## Compilation

If a Makefile is provided, run:

```bash
make
```

Ensure the environment variables for CUDA 11.8 are set correctly:

```bash
export CUDA_HOME=/usr/local/cuda-11.8
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
```

You may need to modify include/library paths in the Makefile depending on your system.

---

## Usage

After compilation, run the executables with properly prepared `.key` files containing all system, bath, and simulation settings.

Example usage (standalone):

```bash
./GPU_2DES key_-100_0.key > out_-100_0.out
```

For automated PBS job submission, refer to:

- `../Run_GPU_2DES.sh` — full batch simulation
- `../Run_GPU_dynamics.sh` — dynamics-only simulations

---

## Notes

- Significantly faster than CPU implementation
- Fully compatible with the ETOM input/output structure and high-level automation scripts
- Designed to run on a **single GPU** per job (multi-GPU not supported in current version)

📄 For the full workflow, refer to the main project README or the following step-by-step documentation:

- `../etom/README.md` — for generating ETOM bath parameters
- `../script/IO-scripts/README.md` — for editing input templates
- `../gen-2d-data/README.md` — for post-processing and visualization
