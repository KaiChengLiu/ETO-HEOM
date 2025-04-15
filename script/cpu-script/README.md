# Job Submission Scripts for 2DES Simulations (PBS + Parallel Execution)

This directory contains two bash scripts to manage and run 2DES simulations using the `CPU_2DES` executable in a PBS-managed HPC environment.

---

## Script: `submit_jobs.sh

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
    - Loops through all combinations of `TAU`, `T`
    - Runs simulations in parallel (up to 11 at a time
    - Outputs standard logs and errors to the `pbslog` and `pbserr` folders
  - Automatically submits the PBS job using `qsub`

### Usage

```bash
./submit_jobs.sh
```

This will create and submit 8 PBS scripts:
- `pbs-script/JOBNAME_1.pbs` to `pbs-script/JOBNAME_8.pbs`
- Logs will be stored in `pbs-script/pbserr/` and `pbs-script/pbslog/`

---

## Summary

| Script            | Purpose                                     | Key Feature                          |
|-------------------|---------------------------------------------|---------------------------------------|
| `submit_jobs.sh`  | Automatically split and submit jobs in bulk | Scales up to multiple TAU chunks     |

These scripts ensure efficient job distribution and error handling in a PBS-based HPC environment.

