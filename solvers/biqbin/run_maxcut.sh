#!/bin/bash
#SBATCH --partition=westmere
#SBATCH --nodes=2
#SBATCH --job-name=biqbin
#SBATCH --output=output_file.txt
#SBATCH --error=error_file.txt

unset LD_PRELOAD
source setupenv.sh

mpirun -n 24 biqbin ./matlab/data/$1 params $2 5 30

