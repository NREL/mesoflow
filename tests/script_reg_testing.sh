#!/bin/bash
TOPDIR=${PWD}

#number of processors to use
NPROCS=${1:-16}

#need mpirun command (on eagle it is srun and 
#peregrine it is mpirun)
MPI_RUN_COMMAND=${2:-srun}

#script that is used to set environment variables for
#post processing or doing conda related commands
SCRIPTFILE=${3:-~/ytenv.sh}

RESULTS_DIR=${4:-~/mflo_results}

MAKEFILE_OPTIONS=${5:-}
rm -rf ${RESULTS_DIR}

declare -a allcases=('Cartesian/Channel' 'Cartesian/ShockTube_N2' 'Cartesian/ShockTube_N2He' 'Cartesian/SpeciesEqTests')

#clean directories
for case in "${allcases[@]}";
do
	cd ${case}
        make realclean
        rm -rf finalpl* *.png
        . ../clean.sh
        cd ${TOPDIR}
done

#run cases
for case in "${allcases[@]}";
do
        echo ${case}
	cd ${case}
        make -j ${MAKEFILE_OPTIONS} 
        . ./runcase.sh "${MPI_RUN_COMMAND}" ${NPROCS} "mflo.order_hyp=5"
        cd ${TOPDIR}
done

source ${SCRIPTFILE}
mkdir ${RESULTS_DIR}

#post process
for case in "${allcases[@]}";
do
	cd ${case}
        . ./verifycase.sh
        mv *.png ${RESULTS_DIR}
        cd ${TOPDIR}
done

cd ${RESULTS_DIR}
convert *.png regtest_results.pdf
