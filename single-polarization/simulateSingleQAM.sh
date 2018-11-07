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


for coreInd in {1..16}
do
RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateSingleQAM($coreInd,32e9,2);\"" & 

sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/run
cp -rf $TMPDIR/* $SLURM_SUBMIT_DIR/run
rm -rf $TMPDIR/*

#End of script

