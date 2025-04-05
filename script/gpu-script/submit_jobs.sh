#!/bin/bash

# Set HOME path
HOME_PATH="YOUR_HOME_PATH"

# Job name (auto set by setup script)
JOBNAME="YOUR_JOBNAME"

# GPU executable
GPU_2DES="$HOME_PATH/etoheom-gpu/src/GPU_2DES"

# Directories
INPUT_DIR="./2d-input"
OUTPUT_DIR="./2d-output"
ERR_DIR="./pbs-script/pbserr"
LOG_DIR="./pbs-script/pbslog"

# TAU and T settings
TLIST="0"
START_TAU=-600
END_TAU=600
STEP_TAU=10

# Create output/error dirs if not exist
mkdir -p "$ERR_DIR"
mkdir -p "$LOG_DIR"

# Generate PBS script
SCRIPT_FILE="$./pbs-script/{JOBNAME}_GPU.pbs"
cat << EOF > "$SCRIPT_FILE"
#!/bin/sh
#PBS -N ${JOBNAME}_GPU
#PBS -l walltime=168:00:00
#PBS -l mem=3GB,nodes=gpu02,ncpus=1
#PBS -q GPU_queue
#PBS -e ${ERR_DIR}/${JOBNAME}_GPU.pbserr
#PBS -o ${LOG_DIR}/${JOBNAME}_GPU.pbslog

source /usr/local/pbs-tools/pbs_prologue.sh

export CUDA_HOME=/usr/local/cuda-11.8
export PATH=\$CUDA_HOME/bin:\$PATH
export LD_LIBRARY_PATH=\$CUDA_HOME/lib64:\$LD_LIBRARY_PATH

cd \$PBS_O_WORKDIR

proc_count=0

for TAU in \$(seq $START_TAU $STEP_TAU $END_TAU); do
    for T in $TLIST; do
        INPUT_FILE="${INPUT_DIR}/key_\${TAU}_\${T}.key"
        OUTPUT_FILE="${OUTPUT_DIR}/out_\${TAU}_\${T}.out"

        if [ ! -f "\$INPUT_FILE" ]; then
            echo "Input file \$INPUT_FILE does not exist!" >> \$PBS_O_WORKDIR/${ERR_DIR}/missing_inputs_GPU.log
            continue
        fi

        ${GPU_2DES} "\$INPUT_FILE" > "\$OUTPUT_FILE" &

        proc_count=\$((proc_count + 1))
        if [ "\$proc_count" -ge 11 ]; then
            wait
            proc_count=0
        fi
    done
done

wait
echo "=== GPU_2DES completed at \$(date) ==="
source /usr/local/pbs-tools/pbs_epilogue.sh
EOF

# Submit the job
qsub "$SCRIPT_FILE"

