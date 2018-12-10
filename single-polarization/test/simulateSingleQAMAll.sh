#!/bin/bash

for i in {1..209}
do

echo "$i"
sbatch --export=RunNo=$1,ParamNo=$i   jobScriptCRLB.txt

done
