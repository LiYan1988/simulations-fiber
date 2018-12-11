#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateScenario203
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH -o simulateScenario203.stdout
#SBATCH -e simulateScenario203.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "-12;0;105000000000;100000000000" "-11;0;105000000000;100000000000" "-10;0;105000000000;100000000000" "-9;0;105000000000;100000000000" "-8;0;105000000000;100000000000" "-7;0;105000000000;100000000000" "-6;0;105000000000;100000000000" "-5;0;105000000000;100000000000" "-4;0;105000000000;100000000000" "-3;0;105000000000;100000000000" "-2;0;105000000000;100000000000" "-1;0;105000000000;100000000000" "0;0;105000000000;100000000000" "1;0;105000000000;100000000000" "2;0;105000000000;100000000000" "3;0;105000000000;100000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario203
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/simulateScenario203
rm -rf $TMPDIR/*

#End of script

