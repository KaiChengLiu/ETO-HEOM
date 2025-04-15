# GPU-Based Job Scripts for 2DES and Population Dynamics Simulations

This directory contains PBS-compatible scripts for running 2DES and population dynamics simulations using GPU acceleration. These scripts are designed for use on a high-performance computing (HPC) cluster with PBS job scheduling and CUDA-enabled GPU nodes.

---

## Script: `Run_GPU_dynamics.sh`

This script runs a population dynamics simulation using the `GPU_dynamics` executable on a GPU node.

### PBS Settings

```bash
#PBS -N GPU_dynamics_job
#PBS -l walltime=168:00:00
#PBS -l mem=3GB,nodes=gpu02,ncpus=1
#PBS -q GPU_queue
#PBS -o Run_GPU_dynamics.pbslog
#PBS -e Run_GPU_dynamics.pbserr
```

### Script Behavior

- Loads CUDA 11.8 environment
- Changes to working directory (`$PBS_O_WORKDIR`)
- Runs the simulation using input: `../src/population.key`
- Output is written to: `../src/dynamics.out`
- Uses prologue/epilogue hooks for job tracking

---

## Script: `submit_jobs.sh`

This is a **job generation and submission script** that automates the creation of a PBS job for running `GPU_2DES` across all `TAU` and `T` combinations.

### Customizable Parameters

```bash
HOME_PATH="/path/to/ETO-HEOM/"                     # Root directory
JOBNAME="JOBNAME"                          # Job name tag
GPU_2DES="$HOME_PATH/etoheom-gpu/src/GPU_2DES"  # Path to executable
INPUT_DIR="./2d-input"                            # Input folder
OUTPUT_DIR="./2d-output"                          # Output folder
ERR_DIR="./pbs-script/pbserr"                     # Error logs
LOG_DIR="./pbs-script/pbslog"                     # Output logs

TLIST="0"               # Population time list
START_TAU=-600          # Starting delay time
END_TAU=600             # Ending delay time
STEP_TAU=10             # Step size for TAU
```

### Cuda Environment Setup

To run GPU jobs, it is necessary to setup cuda environment to access GPU correctly
```bash
export CUDA_HOME=/usr/local/cuda-XX.X
export PATH=\$CUDA_HOME/bin:\$PATH
export LD_LIBRARY_PATH=\$CUDA_HOME/lib64:\$LD_LIBRARY_PATH
```

### Behavior

- Creates missing log directories if needed
- Dynamically generates a PBS script with loops over `TAU` and `T`
- Runs `GPU_2DES`
- Saves PBS error and output logs to `pbserr/` and `pbslog/`
- Submits the job via `qsub`

### Usage

```bash
./submit_jobs.sh
```

Generated script:  
- `${JOBNAME}_GPU.pbs` in `pbs-script/`
- Automatically submitted to the queue

---

## Summary

| Script Name             | Purpose                                    | Notes                                       |
|--------------------------|--------------------------------------------|---------------------------------------------|
| `Run_GPU_dynamics.sh`    | Run population dynamics simulation         | Uses single input key and outputs dynamics  |
| `submit_jobs.sh`         | Auto-generate and submit PBS 2DES job      | Flexible and customizable                   |

These scripts provide a convenient and efficient way to run large-scale GPU-accelerated quantum dynamics simulations on PBS-based HPC systems.

