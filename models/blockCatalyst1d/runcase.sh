#!/bin/bash
# USAGE example: sh runcase.sh mpiexec 8
# run in x direction
$1 -n $2 ./*.ex inputs_x
finfile=$(ls -d plt?????/ | tail -n 1)
rm -rf finalplt_x #Remove existing finalplt_x directory
mv ${finfile} finalplt_x #Move the last plt file for plotting
. ./clean.sh #Remove plt/chk files

# Repeat for other two directions
$1 -n $2 ./*.ex inputs_y
finfile=$(ls -d plt?????/ | tail -n 1)
rm -rf finalplt_y
mv ${finfile} finalplt_y
. ./clean.sh
$1 -n $2 ./*.ex inputs_z
finfile=$(ls -d plt?????/ | tail -n 1)
rm -rf finalplt_z
mv ${finfile} finalplt_z
. ./clean.sh
