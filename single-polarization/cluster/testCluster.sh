#!/bin/bash 
#SBATCH -n 1  

#SBATCH -p glenn 
#SBATCH -A C3SE407-15-3
#SBATCH -o testCluster.stdout
#SBATCH -o testCluster.stderr
#SBATCH -t 0-00:10:00 

cp -pr * $TMPDIR/
cd $TMPDIR

module load MATLAB
matlab -nodesktop -nosplash -singleCompThread -r "testCluster.m"

wait

cp -pr * $SLURM_SUBMIT_DIR/ 
