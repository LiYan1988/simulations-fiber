#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateSingleQAM1
#SBATCH -N 1
#SBATCH -t 0-02:00:00
#SBATCH -o simulateSingleQAM1.stdout
#SBATCH -e simulateSingleQAM1.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "-20;32000000000" "-20;64000000000" "-19;32000000000" "-19;64000000000" "-18;32000000000" "-18;64000000000" "-17;32000000000" "-17;64000000000" "-16;32000000000" "-16;64000000000" "-15;32000000000" "-15;64000000000" "-14;32000000000" "-14;64000000000" "-13;32000000000" "-13;64000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateSingleQAM(${arr[0]},${arr[1]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateSingleQAM1
cp -rf $TMPDIR/* $SLURM_SUBMIT_DIR/simulateSingleQAM1
rm -rf $TMPDIR/*

#End of script

