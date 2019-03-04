#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A C3SE407-15-3
#SBATCH -J simulateScenario541
#SBATCH -n 4
#SBATCH -t 0-10:00:00
#SBATCH -o simulateScenario541.stdout
#SBATCH -e simulateScenario541.stderr


module load MATLAB

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "-7;0;18000000000;50000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario541
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/simulateScenario541
rm -rf $TMPDIR/*

#End of script

