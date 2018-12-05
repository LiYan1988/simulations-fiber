#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateScenario98
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH -o simulateScenario98.stdout
#SBATCH -e simulateScenario98.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "-19;0;26000000000;200000000000" "-18;0;26000000000;200000000000" "-17;0;26000000000;200000000000" "-16;0;26000000000;200000000000" "-15;0;26000000000;200000000000" "-14;0;26000000000;200000000000" "-13;0;26000000000;200000000000" "-12;0;26000000000;200000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario98
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/simulateScenario98
rm -rf $TMPDIR/*

#End of script

