#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J singleQAM
#SBATCH -N 1
#SBATCH -t 0-00:10:00
#SBATCH -o singleQAM.stdout
#SBATCH -e singleQAM.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=( "one;1" "two;2" "three;3" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateSingleQAM(${arr[0]},${arr[1]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/run
cp -rf $TMPDIR/* $SLURM_SUBMIT_DIR/run
rm -rf $TMPDIR/*

#End of script

