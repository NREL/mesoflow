#!/bin/bash
#SBATCH -J yo         # job name
#SBATCH -o out.o%j         # output and error file name (%j expands to jobID)
#SBATCH -N 1               # Total # of nodes
#SBATCH -n 50              # total number of mpi tasks requested
#SBATCH -p development    # queue (partition) -- normal, development, etc.
#SBATCH -t 01:59:00        # run time (hh:mm:ss) - 4.0 hours
#SBATCH --mail-user=ryou@utexas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

#export OMP_NUM_THREADS=48
#./bem_advdiff_adi<advdiff.inp ./out.log  # run the MPI executable named a.out
ibrun ./mesoflow3d.gnu.MPI.ex<inputs ./out.log