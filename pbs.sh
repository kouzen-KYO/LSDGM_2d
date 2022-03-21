#!/bin/sh

# PBS script example for parallel jobs using Open MPI/GNU 4.1
# rkawamot

# sets the number of nodes & processors per node (ppn)
# can also explicitly specify which nodes to use: nodes001:ppn=6+n010 [...]
#PBS -l nodes=001:ppn=12

# name of your job
#PBS -N m_shear_0.2

# send e-mail notifications to these addresses
#PBS -M rkawamot@caltech.edu

# Combine output and error streams into single file
#PBS -j oe

# Name of queue to submit job to
#PBS -q default

# sets the maximum amount of time the job can run for (hr:min:sec)
#PBS -l walltime=24:00:00
numproc=$(cat $PBS_NODEFILE | wc -l)
cd $PBS_O_WORKDIR
echo The number of processes is...$numproc
echo Working directory is...$(pwd)

# modify to include the correct mpirun and your executable
./Main 0.2