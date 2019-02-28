#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A C3SE407-15-3
#SBATCH -J simulateScenario12
#SBATCH -n 3
#SBATCH -t 0-10:00:00
#SBATCH -o simulateScenario12.stdout
#SBATCH -e simulateScenario12.stderr


module load MATLAB

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "-9;0;1000000000;150000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario12
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/simulateScenario12
rm -rf $TMPDIR/*

#End of script

