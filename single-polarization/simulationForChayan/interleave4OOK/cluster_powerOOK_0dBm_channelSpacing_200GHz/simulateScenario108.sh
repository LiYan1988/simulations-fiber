#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateScenario108
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH -o simulateScenario108.stdout
#SBATCH -e simulateScenario108.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "-1;0;28000000000;200000000000" "0;0;28000000000;200000000000" "1;0;28000000000;200000000000" "2;0;28000000000;200000000000" "3;0;28000000000;200000000000" "4;0;28000000000;200000000000" "5;0;28000000000;200000000000" "6;0;28000000000;200000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario108
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/simulateScenario108
rm -rf $TMPDIR/*

#End of script

