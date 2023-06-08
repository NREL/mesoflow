#!/bin/bash
$1 -n $2 ./*.ex inputs_x $3
finfile=$(ls -d plt????? | tail -n 1)
mv ${finfile} finalplt_x
. ../clean.sh
$1 -n $2 ./*.ex inputs_y $3
finfile=$(ls -d plt????? | tail -n 1)
mv ${finfile} finalplt_y
. ../clean.sh
$1 -n $2 ./*.ex inputs_z $3
finfile=$(ls -d plt????? | tail -n 1)
mv ${finfile} finalplt_z
. ../clean.sh
