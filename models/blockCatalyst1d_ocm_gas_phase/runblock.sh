#!/bin/bash

#SBATCH --partition=chem
#SBATCH --time=1:00:00 #Wall Time Limit
#SBATCH --ntasks-per-node=32 #CPU Count
#SBATCH --qos=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lch224@lehigh.edu
#SBATCH --job-name="testMesoflow"
#SBATCH --output="Job%j.%N.out" # Write Standard Output and Error


cd $SLURM_SUBMIT_DIR

# To run jobs, add the following line to your submit script to load modules that are optimized for the underlying CPU (for debug, enge, chem, im2080, health, hawk and infolab partitions)
# source /share/Apps/compilers/etc/lmod/zlmod.sh
module load gcc/9.3.0
module load mvapich2/2.3.4
module load amrex/23.02


# You can change the name for the output file ("output.out")
echo $AMREX_HOME
srun ./mesoflow3d.gnu.MPI.ex inputs_x



exit