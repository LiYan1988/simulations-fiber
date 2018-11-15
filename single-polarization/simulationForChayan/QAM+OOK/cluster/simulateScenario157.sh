#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateScenario157
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH -o simulateScenario157.stdout
#SBATCH -e simulateScenario157.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "8;3;64000000000;100000000000" "9;3;64000000000;100000000000" "10;3;64000000000;100000000000" "-10;4;64000000000;100000000000" "-9;4;64000000000;100000000000" "-8;4;64000000000;100000000000" "-7;4;64000000000;100000000000" "-6;4;64000000000;100000000000" "-5;4;64000000000;100000000000" "-4;4;64000000000;100000000000" "-3;4;64000000000;100000000000" "-2;4;64000000000;100000000000" "-1;4;64000000000;100000000000" "0;4;64000000000;100000000000" "1;4;64000000000;100000000000" "2;4;64000000000;100000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario157
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/simulateScenario157
rm -rf $TMPDIR/*

#End of script

