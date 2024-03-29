#!/bin/bash -l

#SBATCH -D ./
#SBATCH --export=ALL
#SBATCH -p course
#SBATCH -t 5
#SBATCH --cpus-per-task=32
# ADD INSTRUCTIONS TO LOAD THE MODULES HERE
# Example: module load mpi
module load compilers/intel/2019u5 
module load mpi/intel-mpi/2019u5/bin

# ADD COMPILER INSTRUCTION HERE.
# Example: mpicc -o my_mpi_program my_mpi_program.c
mpicc -std=c99 -fopenmp main-mpi.c coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c -o mpi -lm


# SLURM_NTASKS is given by the -n flag when using sbatch. 
procs=${SLURM_NTASKS:-1}
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#echo using ${OMP_NUM_THREADS} OpenMP threads
# RUN THE PROGRAM HERE USING $procs WITH THE -np FLAG.
# Example: mpirun -np $procs ./my_mpi_program
mpirun -np $procs ./mpi 512_coords.coord cout fout nout
