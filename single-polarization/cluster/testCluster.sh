#!/bin/bash 
#SBATCH -N 1 -n 1  

#SBATCH -p glenn 
#SBATCH -A C3SE2018-1-15
#SBATCH -o testCluster.stdout
#SBATCH -e testCluster.stderr
#SBATCH -t 0-00:10:00 

echo $TMPDIR
echo $SLURM_SUBMIT_DIR

cp -pr * $TMPDIR/
cd $TMPDIR

module load matlab
matlab -nodesktop -nosplash -singleCompThread -r "testCluster"

wait

cp -pr * $SLURM_SUBMIT_DIR/ 

echo $TMPDIR
echo $SLURM_SUBMIT_DIR