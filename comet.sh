#!/bin/bash

#SBATCH -A TG-ENG160025
#SBATCH --job-name="cyc_compress"
#SBATCH --output="cyc_compress.%j"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -t 06:00:00
#SBATCH --export=ALL

export OMP_NUM_THREADS=12
ibrun -np 1 ./MainIsoBox $1 $2 