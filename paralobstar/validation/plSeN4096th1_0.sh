#!/bin/bash
#PBS -N plSeN4096th1_0
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:20:00
#PBS -l pmem=1gb
#PBS -q tiny
#PBS -m aeb -M johannes-stefan.martin@student.uni-tuebingen.de

source ~/.bashrc

export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH

# Loading modules
module load mpi/openmpi/3.1-gnu-9.2
module load lib/hdf5/1.10.7-openmpi-3.1-gnu-9.2

# Going to working directory
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

# creating output directory
mkdir -p $PBS_JOBNAME

# Starting program
mpirun --bind-to core --map-by core -report-bindings ../bin/paralobstar -c config/$PBS_JOBNAME.info -p $PBS_JOBNAME.h5



