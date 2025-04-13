# INSTALLATION & JOB SETUP GUIDE

This document provides system requirements, build instructions, and job setup procedures for running 2DES simulations using the ETO-HEOM framework.

---

## System Requirements

The code is developed for Linux systems with support for both CPU and NVIDIA GPU execution.

### Required Software

| Component      | Version    |
|----------------|------------|
| C++ Standard   | C++14      |
| GSL            | 2.8        |
| CUDA           | 11.8       |
| Python         | 3.11.4     |
| NumPy          | 1.24.3     |
| SciPy          | 1.11.1     |
| Matplotlib     | 3.7.1      |


> GSL is required to CPU HEOM solvers calculation.
> CUDA is required to GPU HEOM solvers calculation.
> Python and associated libraries are required to generate visualize results.

---

## Build Instructions

Before setting up any jobs, first check the compiler configuration to fit your system,
you can modify these variables accordingly:
```make
CXX := g++
CXXFLAGS := -Wall -Wno-sign-compare -Wno-endif-labels -Wno-unused-variable -O3 -std=c++14

# CUDA compiler and flags
NVCC := $/path/to/cuda-XX.X/bin/nvcc
NVCCFLAGS := -g -Xcompiler="-Wall -Wno-sign-compare -Wno-endif-labels -Wno-unused-variable" -O3 -std=c++14

# Include and library paths
INCLUDES := -I/usr/include
CUDA_INC := -I$/path/to/cuda-XX.X/include
CUDA_LIB := -L$/path/to/cuda-XX.X/lib64
```

After compilation, the following binaries will be generated in their respective directories:
```bash
etoheom-cpu/src/CPU_2DES
etoheom-cpu/src/CPU_dynamics
etoheom-gpu/src/GPU_2DES
etoheom-gpu/src/GPU_dynamics
```

Now you can compile the ETO-HEOM using the unified Makefile in the project root:
```bash
cd /path/to/ETO-HEOM/
make
```

If you only want to build one version:
```bash
make cpu        # Build only CPU version
make gpu        # Build only GPU version
```

To clean all build artifacts:
```bash
make clean      # Clean all object files and binaries
```
---

## Job Setup Scripts

### CPU Job Setup

Initializes a simulation folder for CPU-based 2DES simulations.

**Usage:**
```bash
./setup_cpu_job.sh JOBNAME BATHTYPE
```

**Arguments:**
- `JOBNAME`: Custom identifier for the simulation (e.g., `dimer`)
- `BATHTYPE`: Type of spectral density — must be one of:
  - `debye_lorentz`
  - `ohmic`
  - `superohmic`

**Resulting Structure:**
```
CPU_JOBNAME_BATHTYPE/
├── 2d-input/                            # Input .key files
├── 2d-output/                           # Output .out files
├── pbs-script/                          # PBS scripts
│   ├── pbserr/                          # PBS error logs
│   └── pbslog/                          # PBS output logs
├── key.key-tmpl                         # Input template file
├── gen_2d_data.sh                       # Script to generate .key files
├── clean_2d_data.sh                     # Script to clean .key and .out files
├── gen_2d_spectrum.py                   # Script to gen 2d spectrum
├── gen_2d_spectrum_example.ipynb        # # gen 2d spectrum example
├── BATHTYPE_ETOM.py                     # ETO model for the specified bath type
├── BATHTYPE_ETOM_example.ipynb          # ETO model example for the specified bath type
├── submit_jobs.sh                       # Script to generate & submit PBS jobs (CPU)
└── README.md                            # Usage of working file

```

---

### GPU Job Setup

Initializes a simulation folder for GPU-based 2DES simulations.

**Usage:**
```bash
./setup_gpu_job.sh JOBNAME BATHTYPE
```

**Arguments:**
- `JOBNAME`: Custom identifier for the simulation (e.g., `dimer`)
- `BATHTYPE`: Type of spectral density — must be one of:
  - `debye_lorentz`
  - `ohmic`
  - `superohmic`

**Resulting Structure:**
```
GPU_JOBNAME_BATHTYPE/
├── 2d-input/                            # Input .key files
├── 2d-output/                           # Output .out files
├── pbs-script/                          # PBS scripts
│   ├── pbserr/                          # PBS error logs
│   └── pbslog/                          # PBS output logs
├── key.key-tmpl                         # Input template file
├── gen_2d_data.sh                       # Script to generate .key files
├── clean_2d_data.sh                     # Script to clean .key and .out files
├── gen_2d_spectrum.py                   # Script to analyze output
├── BATHTYPE_ETOM.py                     # ETO model for the specified bath type
├── BATHTYPE_ETOM_example.ipynb          # ETO model example for the specified bath type
├── submit_jobs.sh                       # Script to generate & submit PBS jobs (GPU)
└── README.md                            # Usage of working file
```

---

## Data Generation and Cleanup

### Generate `.key` Files

After setup, run the following to generate input files:
```bash
./gen_2d_data.sh
```

This script creates `.key` files inside the `2d-input/` directory, covering all combinations of **TAU** and **T** values.

---

### Clean `.key` and `.out` Files

To remove all generated input and output files:
```bash
./clean_2d_data.sh
```

This deletes:

- `2d-input/*.key`
- `2d-output/*.out` 
