unset LD_PRELOAD
module load OpenBLAS 
module load OpenMPI/4.1.4-GCC-11.3.0
module load imkl
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
