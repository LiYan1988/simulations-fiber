#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateScenario91
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH -o simulateScenario91.stdout
#SBATCH -e simulateScenario91.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "2;-5;64000000000;50000000000" "3;-5;64000000000;50000000000" "4;-5;64000000000;50000000000" "5;-5;64000000000;50000000000" "6;-5;64000000000;50000000000" "7;-5;64000000000;50000000000" "8;-5;64000000000;50000000000" "9;-5;64000000000;50000000000" "10;-5;64000000000;50000000000" "-10;-4;64000000000;50000000000" "-9;-4;64000000000;50000000000" "-8;-4;64000000000;50000000000" "-7;-4;64000000000;50000000000" "-6;-4;64000000000;50000000000" "-5;-4;64000000000;50000000000" "-4;-4;64000000000;50000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario91
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/simulateScenario91
rm -rf $TMPDIR/*

#End of script

