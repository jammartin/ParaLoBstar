#!/bin/bash
#PBS -N plummerN1e4Hnp10
#PBS -l nodes=1:ppn=10
#PBS -l walltime=00:05:00
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

# Starting program
mpirun --bind-to core --map-by core -report-bindings bin/paralobstar



