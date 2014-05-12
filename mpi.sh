#!/bin/sh
#SBATCH -p short
#SBATCH -n16

R_LIBS=/scratch/lustre/mare9467/lib64/R/library
export R_LIBS

mpirun --vanilla --slave -f /scratch/lustre/home/mare9467/mif_home/RRR/one.R
