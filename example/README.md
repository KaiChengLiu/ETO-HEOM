# DIMER Example (CPU & GPU Versions)

This directory contains example working files for simulating a **dimer system** using the ETO-HEOM framework, including both **CPU** and **GPU** versions.

---

## Overview

The dimer example is designed to help you understand the complete workflow of using **ETO-HEOM** to generate a **2D electronic spectrum**. It includes all necessary files and scripts to simulate and visualize the response of a dimer system in a structured bath environment.

---

## Contents

- `CPU_dimer_superohmic/`  
  Example files configured to run with the CPU version of ETO-HEOM (`CPU_2DES`).

- `GPU_dimer_superohmic/`  
  Example files configured to run with the GPU version of ETO-HEOM (`GPU_2DES`).

Each folder contains:
- Input template (`key.key-tmpl`)
- ETOM bath generation script (`{BATHTYPE}_ETOM.py`)
- Data generation script (`gen_2d_data.sh`)
- Submission scripts (`submit_jobs.sh` or `Run_GPU_2DES.sh`)
- Visualization script (`gen_2d_data.py`)
- Example output or plotting files

---

## How to Use

Follow the step-by-step instructions provided in each folder’s `README.md` to:

1. **Generate ETOM parameters** based on the bath type
2. **Set up input files** for the HEOM simulation
3. **Run the simulations** using either CPU or GPU
4. **Visualize the resulting 2D spectrum**

These examples serve as a template for your own systems. You can modify the Hamiltonian, bath model, dipole orientations, and simulation parameters as needed.

---

## Recommended Starting Point

To begin practicing with ETO-HEOM, refer to the `README.md` inside:

- `./cpu/` for CPU version
- `./gpu/` for GPU version

Each README guides you through the simulation pipeline from input setup to final spectrum generation.

---

## Requirements

- For **CPU version**: GSL 2.8 and a C++ compiler
- For **GPU version**: CUDA 11.8 and compatible NVIDIA GPU

---

These examples are a great entry point for users who want to learn how to apply ETO-HEOM to real systems and generate high-resolution 2D spectra.

