# Job Submission Scripts for 2DES Simulations (PBS + Parallel Execution)

This directory contains two bash scripts to manage and run 2DES simulations using the `CPU_2DES` executable in a PBS-managed HPC environment.

---

## Script 1: `Run_CPU_2DES.sh`

This script is designed to **run a batch of 2DES simulations** on a single node using **up to 11 parallel processes**.

### PBS Directives

```bash
#PBS -N Run_CPU_2DES         # Job name
#PBS -r n                    # Do not rerun job if it fails
#PBS -e Run_CPU_2DES.pbserr  # File for standard error output
#PBS -o Run_CPU_2DES.pbslog  # File for standard output log
#PBS -q long                 # Queue name
#PBS -l walltime=720:00:00,nodes=1:ppn=11  # Resource request: 1 node with 11 processors
```

### Parameters and Logic

- `CPU_2DES`: path to the simulation executable
- `TLIST`: list of population times T
- `TAU`: ranges from -600 to 600 with step size 10
- For each combination of `TAU` and `T`, it checks if the input file exists (`key_{TAU}_{T}_wc125.key`), and runs the simulation if so.
- Executes up to 11 parallel jobs at a time using `&` and `wait`.

### Output

- Results are written to `../2d_output/out_{TAU}_{T}_wc125.out`
- Errors from missing inputs are appended to `error.log`
- PBS logs and errors go to `Run_CPU_2DES.pbslog` and `Run_CPU_2DES.pbserr`

---

## Script 2: `submit_jobs.sh`

This is a **meta-script** that automatically generates and submits **multiple PBS scripts** to run 2DES simulations in parallel batches across different `TAU` values.

### Parameters

```bash
HOME_PATH="/home/andrew91411"         # User's home path
JOBNAME="$1"                          # Name prefix for each job (passed as script argument)
CPU_2DES="$HOME_PATH/ETO-HEOM/etoheom-cpu/src/CPU_2DES"  # Path to the simulation executable

INPUT_DIR="./2d-input"               # Input files directory
OUTPUT_DIR="./2d-output"             # Output files directory
PBS_DIR="."                          # PBS scripts will be placed here
ERR_DIR="./pbs-script/pbserr"        # Error log directory
LOG_DIR="./pbs-script/pbslog"        # Output log directory
```

```bash
TLIST="0"         # Population time list
FLIST="100"       # Dummy loop variable (e.g., pulse width or field strength)
START_TAU=-600    # Start of TAU scan
END_TAU=600       # End of TAU scan
STEP_TAU=10       # Step of TAU
NUM_SCRIPTS=8     # Number of PBS jobs to submit (splits TAU range into 8 parts)
```

### How it Works

- Computes how many `TAU` values each script should handle
- For each batch:
  - Generates a PBS script (`${JOBNAME}_${i}.pbs`) with the correct TAU range
  - Each script:
    - Loops through all combinations of `TAU`, `T`, and `F`
    - Runs simulations in parallel (up to 11 at a time)
    - Handles missing input files by logging to `missing_inputs_{i}.log`
    - Outputs standard logs and errors to the `pbslog` and `pbserr` folders
  - Automatically submits the PBS job using `qsub`

### Usage

```bash
bash submit_jobs.sh <JOBNAME>
```

Example:

```bash
bash submit_jobs.sh my2desrun
```

This will create and submit 8 PBS scripts:
- `pbs-script/my2desrun_1.pbs` to `pbs-script/my2desrun_8.pbs`
- Logs will be stored in `pbs-script/pbserr/` and `pbs-script/pbslog/`

---

## Summary

| Script            | Purpose                                     | Key Feature                          |
|-------------------|---------------------------------------------|---------------------------------------|
| `Run_CPU_2DES.sh` | Run 2DES jobs on a single node              | Manages 11 parallel processes         |
| `submit_jobs.sh`  | Automatically split and submit jobs in bulk | Scales up to multiple TAU chunks     |

These scripts ensure efficient job distribution and error handling in a PBS-based HPC environment.

