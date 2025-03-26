#!/bin/bash

# Set the path to the CPU_2DES program
CPU_2DES="/home/andrew91411/ETO-HEOM/etoheom-cpu/src/CPU_2DES"
PBS_DIR="."  
OUTPUT_DIR="./2d_output"
INPUT_DIR="./2d_input"
TLIST="0"
FLIST="100"
START_TAU=-600
END_TAU=600
STEP_TAU=10
NUM_SCRIPTS=8  

# Calculate the number of TAU values per script
TOTAL_TAU=$(( (END_TAU - START_TAU) / STEP_TAU + 1 ))
TAU_PER_SCRIPT=$(( TOTAL_TAU / NUM_SCRIPTS ))

mkdir -p "$PBS_DIR"

# Generate PBS scripts
for i in $(seq 1 $NUM_SCRIPTS); do
    START=$(( START_TAU + (i - 1) * TAU_PER_SCRIPT * STEP_TAU ))
    END=$(( START + (TAU_PER_SCRIPT - 1) * STEP_TAU ))
    if [ "$i" -eq "$NUM_SCRIPTS" ]; then
        END=$END_TAU  
    fi

    #echo $START $END

    SCRIPT_FILE="$PBS_DIR/Run_CPU_2DES_part_${i}.pbs"
    cat << EOF > "$SCRIPT_FILE"
#!/bin/sh
#PBS -N Run_CPU_2DES_part_${i}
#PBS -r n
#PBS -e ./pbserr/Run_CPU_2DES_part_${i}.pbserr
#PBS -o ./pbslog/Run_CPU_2DES_part_${i}.pbslog
#PBS -q long
#PBS -l walltime=720:00:00,nodes=1:ppn=11

# Change to the directory where the job was submitted
cd \$PBS_O_WORKDIR
s
# Loop through TAU values for this script
for TAU in \$(seq $START $STEP_TAU $END); do
    for T in $TLIST; do
        for F in $FLIST; do
          INPUT_FILE="${INPUT_DIR}/key_\${TAU}_\${T}.key"
          OUTPUT_FILE="${OUTPUT_DIR}/out_\${TAU}_\${T}.out"
  
          # Check if the input file exists
          if [ ! -f "\$INPUT_FILE" ]; then
              echo "Input file \$INPUT_FILE does not exist!" >> \$PBS_O_WORKDIR/error_part_${i}.log
              continue
          fi
  
          # Execute the program in the background
          ${CPU_2DES} "\$INPUT_FILE" > "\$OUTPUT_FILE" &
  
          # Increment the process count
          proc_count=\$((proc_count + 1))
  
          # Check if we've started 11 processes, if so, wait for them to finish
          if [ "\$proc_count" -ge 11 ]; then
              wait
              proc_count=0  # Reset the process count after waiting
          fi
        done
    done
done

# Wait for any remaining background processes to finish
wait

echo "All tasks in part ${i} are completed."
EOF

    # Submit the script
    qsub "$SCRIPT_FILE"
done

echo "All PBS scripts have been submitted."

