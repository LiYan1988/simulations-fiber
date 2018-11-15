#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateScenario127
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH -o simulateScenario127.stdout
#SBATCH -e simulateScenario127.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "-10;2;32000000000;100000000000" "-9;2;32000000000;100000000000" "-8;2;32000000000;100000000000" "-7;2;32000000000;100000000000" "-6;2;32000000000;100000000000" "-5;2;32000000000;100000000000" "-4;2;32000000000;100000000000" "-3;2;32000000000;100000000000" "-2;2;32000000000;100000000000" "-1;2;32000000000;100000000000" "0;2;32000000000;100000000000" "1;2;32000000000;100000000000" "2;2;32000000000;100000000000" "3;2;32000000000;100000000000" "4;2;32000000000;100000000000" "5;2;32000000000;100000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario127
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/simulateScenario127
rm -rf $TMPDIR/*

#End of script

