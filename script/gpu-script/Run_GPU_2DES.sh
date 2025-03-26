#!/bin/sh
#PBS -N GPU_2DES_job
#PBS -l walltime=168:00:00
#PBS -l mem=3GB,nodes=gpu02,ncpus=1
#PBS -q GPU_queue
#PBS -e Run_GPU_2DES.pbserr
#PBS -o Run_GPU_2DES.pbslog



source /usr/local/pbs-tools/pbs_prologue.sh

# Load CUDA 11.8 environment
export CUDA_HOME=/usr/local/cuda-11.8
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH

# Define the GPU_2DES executable path
GPU_2DES="../src/GPU_2DES"  # The executable is placed in the PBS working directory

# Change to the directory where the job was submitted
cd $PBS_O_WORKDIR

# Define TLIST
TLIST="0"

# Log job start time
echo "=== Starting GPU_2DES at $(date) ==="

# Loop through all TAU and T combinations
for TAU in $(seq -600 10 600); do
    for T in $TLIST; do
        INPUT_FILE="../2d_input/key_${TAU}_${T}.key"
        OUTPUT_FILE="../2d_output/out_${TAU}_${T}.out"

        # Check if the input file exists
        if [ ! -f "$INPUT_FILE" ]; then
            echo "Input file $INPUT_FILE does not exist!" >> $PBS_O_WORKDIR/error.log
            continue
        fi

        # Execute GPU_2DES and save output to the corresponding file
        ${GPU_2DES} ${INPUT_FILE} > ${OUTPUT_FILE} &

        # Track the number of background processes
        proc_count=$((proc_count + 1))

        # If 11 processes are running, wait for them to finish before continuing
        if [ "$proc_count" -ge 11 ]; then
            wait
            proc_count=0  # Reset process count
        fi
    done
done

# Wait for any remaining background processes to finish
wait

# Log job completion time
echo "=== GPU_2DES execution completed at $(date) ==="
source /usr/local/pbs-tools/pbs_epilogue.sh
