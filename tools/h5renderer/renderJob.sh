#!/bin/bash
#PBS -N ParaLoBstar-h5renderer
#PBS -l nodes=1:ppn=20
#PBS -l walltime=00:05:00
#PBS -l mem=4gb
#PBS -q tiny
#PBS -m aeb -M johannes-stefan.martin@student.uni-tuebingen.de

source ~/.bashrc

# Loading modules
module load lib/boost/1.76.0
module load lib/hdf5/1.10.7-gnu-9.2

# Going to working directory
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

# Defining OpenMP number of threads
export OMP_NUM_THREADS=20

# setting memory allocation limit to 200MB per thread
ulimit -s 200000

# executing program
bin/h5renderer
