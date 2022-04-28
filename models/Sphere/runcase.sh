#!/bin/bash
$1 -n $2 ./*.ex inputs
finfile=$(ls -dt plt*/ | head -n 1)
mv ${finfile} finalplt
. ../clean.sh
