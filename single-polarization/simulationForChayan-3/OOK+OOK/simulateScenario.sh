#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A C3SE407-15-3
#SBATCH -J singleQAM
#SBATCH -n 1
#SBATCH -t 0-00:10:00
#SBATCH -o singleQAM.stdout
#SBATCH -e singleQAM.stderr


module load MATLAB

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=( "one;1;a;r" "two;2;b;s" "three;3;c;t" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]} ${arr[4]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]},${arr[4]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/run
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/run
rm -rf $TMPDIR/*

#End of script

