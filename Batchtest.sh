#!/bin/bash -l

#SBATCH -p course -t 30
#SBATCH --cpus-per-task=32


echo "Node list                    : $SLURM_JOB_NODELIST"
echo "Number of nodes allocated    : $SLURM_JOB_NUM_NODES or $SLURM_NNODES"
echo "Number of threads or processes          : $SLURM_NTASKS"
echo "Number of processes per node : $SLURM_TASKS_PER_NODE"
echo "Requested tasks per node     : $SLURM_NTASKS_PER_NODE"
echo "Requested CPUs per task      : $SLURM_CPUS_PER_TASK"
echo "Scheduling priority          : $SLURM_PRIO_PROCESS"

# Defaults on Barkla (but set to be safe)
#SBATCH -D ./
#SBATCH --export=ALL

module load compilers/intel/2019u5 
module load mpi/intel-mpi/2019u5/bin
# SLURM settings echo
# [Echo commands]

EXE="./MSerial2" 


# Set the number of OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

COORDINATE_FILENAME=$1
OUTPUT_FILENAME=$2

echo "Running ${EXE} with input file ${COORDINATE_FILENAME} and output file ${OUTPUT_FILENAME}"
${EXE} ${COORDINATE_FILENAME} ${OUTPUT_FILENAME}

echo using ${OMP_NUM_THREADS} OpenMP threads




#!/bin/bash



# Array of thread counts
#threads=(1 2 4 8 16 32)

# Loop over thread counts
#for t in "${threads[@]}"
#do
    #echo "Running with $t threads..."
    #export OMP_NUM_THREADS=$t
    #./fomp.exe > "output_$t.txt"
#done

#echo "Done."