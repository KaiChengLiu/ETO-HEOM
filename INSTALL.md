# INSTALLATION & JOB SETUP GUIDE

This document provides system requirements, build instructions, and job setup procedures for running 2DES simulations using the ETO-HEOM framework.

---

## System Requirements

The code is developed for Linux systems with support for both CPU and NVIDIA GPU execution.

### Required Software

| Component      | Version    |
|----------------|------------|
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

Before setting up any jobs, compile the solvers by running `make` in the project root:

```bash
cd /path/to/ETO-HEOM/
make
```

This builds the following binaries:

- **CPU_2DES** — CPU-based HEOM 2DES simulator  
- **GPU_2DES** — GPU-accelerated HEOM 2DES simulator (if applicable)

Ensure the binaries are placed in their respective `src/` directories, e.g.:

```
etoheom-cpu/src/
etoheom-gpu/src/
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
JOBNAME_BATHTYPE/
├── 2d-input/                 # Input .key files
├── 2d-output/                # Output .out files
├── pbs-script/               # PBS scripts
│   ├── pbserr/               # PBS error logs
│   └── pbslog/               # PBS output logs
├── key.key-tmpl              # Input template file
├── gen_2d_data.sh            # Script to generate .key files
├── clean_2d_data.sh          # Script to clean .key and .out files
├── gen_2d_spectrum.py        # Script to analyze output
├── BATHTYPE_ETOM.py          # ETO model for the specified bath type
└── submit_jobs.sh            # Script to generate & submit PBS jobs (CPU)

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
JOBNAME_BATHTYPE/
├── 2d-input/                 # Input .key files
├── 2d-output/                # Output .out files
├── pbs-script/               # PBS scripts
│   ├── pbserr/               # PBS error logs
│   └── pbslog/               # PBS output logs
├── key.key-tmpl              # Input template file
├── gen_2d_data.sh            # Script to generate .key files
├── clean_2d_data.sh          # Script to clean .key and .out files
├── gen_2d_spectrum.py        # Script to analyze output
├── BATHTYPE_ETOM.py          # ETO model for the specified bath type
└── submit_jobs.sh            # Script to generate & submit PBS jobs (GPU)
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
