#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateScenario56
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH -o simulateScenario56.stdout
#SBATCH -e simulateScenario56.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "3;0;18000000000;150000000000" "4;0;18000000000;150000000000" "5;0;18000000000;150000000000" "6;0;18000000000;150000000000" "7;0;18000000000;150000000000" "8;0;18000000000;150000000000" "9;0;18000000000;150000000000" "10;0;18000000000;150000000000" "-20;0;19000000000;150000000000" "-19;0;19000000000;150000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario56
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/simulateScenario56
rm -rf $TMPDIR/*

#End of script

