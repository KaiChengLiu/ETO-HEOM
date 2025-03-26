#!/bin/sh
#PBS -N GPU_dynamics_job
#PBS -l walltime=168:00:00
#PBS -l mem=3GB,nodes=gpu02,ncpus=1
#PBS -q GPU_queue
#PBS -o Run_GPU_dynamics.pbslog
#PBS -e Run_GPU_dynamics.pbserr

source /usr/local/pbs-tools/pbs_prologue.sh

export CUDA_HOME=/usr/local/cuda-11.8
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH

cd $PBS_O_WORKDIR

../src/GPU_dynamics ../src/population.key > ../src/dynamics.out


source /usr/local/pbs-tools/pbs_epilogue.sh
