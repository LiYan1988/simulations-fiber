#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateScenario10
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH -o simulateScenario10.stdout
#SBATCH -e simulateScenario10.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "8;0;3000000000;150000000000" "9;0;3000000000;150000000000" "10;0;3000000000;150000000000" "-20;0;4000000000;150000000000" "-19;0;4000000000;150000000000" "-18;0;4000000000;150000000000" "-17;0;4000000000;150000000000" "-16;0;4000000000;150000000000" "-15;0;4000000000;150000000000" "-14;0;4000000000;150000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario10
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/simulateScenario10
rm -rf $TMPDIR/*

#End of script

