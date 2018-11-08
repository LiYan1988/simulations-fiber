#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateSingleQAM3
#SBATCH -N 1
#SBATCH -t 0-02:00:00
#SBATCH -o simulateSingleQAM3.stdout
#SBATCH -e simulateSingleQAM3.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "-4;32000000000" "-4;64000000000" "-3;32000000000" "-3;64000000000" "-2;32000000000" "-2;64000000000" "-1;32000000000" "-1;64000000000" "0;32000000000" "0;64000000000" "1;32000000000" "1;64000000000" "2;32000000000" "2;64000000000" "3;32000000000" "3;64000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateSingleQAM(${arr[0]},${arr[1]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateSingleQAM3
cp -rf $TMPDIR/* $SLURM_SUBMIT_DIR/simulateSingleQAM3
rm -rf $TMPDIR/*

#End of script

