#!/bin/bash

# Set HOME path prefix
HOME_PATH="/home/andrew91411"

# Set the job name (auto-inserted by setup script)
JOBNAME="CPU_dimer_superohmic"

CPU_2DES="$HOME_PATH/ETO-HEOM/etoheom-cpu/src/CPU_2DES"
PBS_DIR="."  # Current directory
OUTPUT_DIR="./2d-output"
INPUT_DIR="./2d-input"
ERR_DIR="./pbs-script/pbserr"
LOG_DIR="./pbs-script/pbslog"

TLIST="0"
FLIST="100"
START_TAU=-600
END_TAU=600
STEP_TAU=10
NUM_SCRIPTS=8

TOTAL_TAU=$(( (END_TAU - START_TAU) / STEP_TAU + 1 ))
TAU_PER_SCRIPT=$(( TOTAL_TAU / NUM_SCRIPTS ))

mkdir -p "$ERR_DIR"
mkdir -p "$LOG_DIR"

for i in $(seq 1 $NUM_SCRIPTS); do
    START=$(( START_TAU + (i - 1) * TAU_PER_SCRIPT * STEP_TAU ))
    END=$(( START + (TAU_PER_SCRIPT - 1) * STEP_TAU ))
    if [ "$i" -eq "$NUM_SCRIPTS" ]; then
        END=$END_TAU
    fi

    SCRIPT_FILE="${PBS_DIR}/pbs-script/${JOBNAME}_${i}.pbs"
    cat << EOF > "$SCRIPT_FILE"
#!/bin/sh
#PBS -N ${JOBNAME}_${i}
#PBS -r n
#PBS -e ${ERR_DIR}/${JOBNAME}_${i}.pbserr
#PBS -o ${LOG_DIR}/${JOBNAME}_${i}.pbslog
#PBS -q long
#PBS -l walltime=720:00:00,nodes=1:ppn=11

cd \$PBS_O_WORKDIR

proc_count=0

for TAU in \$(seq $START $STEP_TAU $END); do
    for T in $TLIST; do
        for F in $FLIST; do
            INPUT_FILE="${INPUT_DIR}/key_\${TAU}_\${T}.key"
            OUTPUT_FILE="${OUTPUT_DIR}/out_\${TAU}_\${T}.out"

            if [ ! -f "\$INPUT_FILE" ]; then
                echo "Input file \$INPUT_FILE does not exist!" >> \$PBS_O_WORKDIR/${ERR_DIR}/missing_inputs_${i}.log
                continue
            fi

            ${CPU_2DES} "\$INPUT_FILE" > "\$OUTPUT_FILE" &

            proc_count=\$((proc_count + 1))
            if [ "\$proc_count" -ge 11 ]; then
                wait
                proc_count=0
            fi
        done
    done
done

wait
echo "All tasks in part ${i} are completed."
EOF

    qsub "$SCRIPT_FILE"
done

echo "All PBS scripts have been submitted."
