#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateSingleQAM2
#SBATCH -N 1
#SBATCH -t 0-02:00:00
#SBATCH -o simulateSingleQAM2.stdout
#SBATCH -e simulateSingleQAM2.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "-12;32000000000" "-12;64000000000" "-11;32000000000" "-11;64000000000" "-10;32000000000" "-10;64000000000" "-9;32000000000" "-9;64000000000" "-8;32000000000" "-8;64000000000" "-7;32000000000" "-7;64000000000" "-6;32000000000" "-6;64000000000" "-5;32000000000" "-5;64000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateSingleQAM(${arr[0]},${arr[1]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateSingleQAM2
cp -rf $TMPDIR/* $SLURM_SUBMIT_DIR/simulateSingleQAM2
rm -rf $TMPDIR/*

#End of script

