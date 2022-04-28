#!/bin/bash
$1 -n $2 ./*.ex inputs
finfile=$(ls -d plt?????/ | tail -n 1)
mv ${finfile} finalplt
. ../clean.sh
