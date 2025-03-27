# ETO-HEOM (CPU Version)

This directory contains the **CPU implementation** of the Extended Thiele Oscillator Hierarchical Equations of Motion (**ETO-HEOM**) framework for simulating open quantum systems, especially applied to 2D electronic spectroscopy.

---

## Structure

- **`lib/`**  
  Contains the static library (`libcpu2des.a`) and compiled objects for CPU-based HEOM propagation.

- **`src/`**  
  Contains the main executable sources for:
  - `CPU_2DES`: for 2D electronic spectroscopy simulations
  - `CPU_dynamics` or similar: for population dynamics or other time-evolution tasks

---

## Dependencies

- This CPU version of ETO-HEOM uses the **GNU Scientific Library (GSL)** version **2.8** for numerical integration and matrix operations.
- Make sure GSL 2.8 is correctly installed and accessible on your system before compiling or running the code.

---

## Compilation

Compilation is handled via Makefile (if provided), and typically builds the library and main executables. You may need to adjust include/library paths for your system’s GSL installation.

Example (if Makefile is present):

```bash
make
```

---

## Usage

Once compiled, the main executables (e.g., `CPU_2DES`) can be called with a `.key` input file that contains system Hamiltonian, bath parameters, dipole information, and time settings.

Refer to higher-level scripts and documentation in:

- `../script/IO-scripts/README.md` — input structure
- `../Run_CPU_2DES.sh` — PBS job submission

---

## Notes

- This implementation is single-GPU/multi-core CPU optimized but **does not require GPU hardware**.
- Useful for benchmarking, debugging, or running HEOM simulations on CPU-based HPC clusters.
