#!/bin/bash

# Set initial variables
t0=50
propagate_time=600

tau_step=10
tau_bound=600
input_file="key.key-tmpl"
T=(0)  # Use a bash array to store values for T

# Read the file content into a variable once
input_content=$(cat "$input_file")

# Function to process file using awk for text substitution and replacement
process_file() {
    local tau1=$1
    local tau2=$2
    local tau3=$3
    local t_end=$4
    local of_name=$5

    echo "Processing file: $of_name with tau1=$tau1, tau2=$tau2, tau3=$tau3, t_end=$t_end"

    # Use awk for replacement
    echo "$input_content" | awk -v tau1="$tau1" -v tau2="$tau2" -v tau3="$tau3" -v t_end="$t_end" '
    {
        gsub(/TAU1/, tau1);
        gsub(/TAU2/, tau2);
        gsub(/TAU3/, tau3);
        gsub(/T_END/, t_end);
        print;
    }' > "$of_name"
}

# First loop: i from 0 to tau_bound
for t in "${T[@]}"; do
    for ((i=0; i<=tau_bound; i+=tau_step)); do
        tau1=$t0
        tau2=$(echo "$t0 + $i" | bc)
        tau3=$(echo "$t0 + $i + $t" | bc)
        t_end=$(echo "$t0 + $i + $t + $propagate_time" | bc)

        of_name="key_${i}_${t}.key"
        process_file "$tau1" "$tau2" "$tau3" "$t_end" "$of_name" &
    done
done

# Second loop: i from -tau_bound to 0
for t in "${T[@]}"; do
    for ((i=-tau_bound; i<0; i+=tau_step)); do
        abs_i=$(echo "${i#-}")  # Get the absolute value of i
        tau1=$(echo "$t0 + $abs_i" | bc)
        tau2=$t0
        tau3=$(echo "$t0 + $abs_i + $t" | bc)
        t_end=$(echo "$t0 + $abs_i + $t + $propagate_time" | bc)

        of_name="key_${i}_${t}.key"
        process_file "$tau1" "$tau2" "$tau3" "$t_end" "$of_name" &
    done
done

# Wait for all background processes to complete
wait
