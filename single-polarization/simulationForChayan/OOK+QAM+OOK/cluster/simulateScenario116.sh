#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateScenario116
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH -o simulateScenario116.stdout
#SBATCH -e simulateScenario116.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "3;-7;32000000000;150000000000" "4;-7;32000000000;150000000000" "5;-7;32000000000;150000000000" "6;-7;32000000000;150000000000" "7;-7;32000000000;150000000000" "8;-7;32000000000;150000000000" "9;-7;32000000000;150000000000" "10;-7;32000000000;150000000000" "-10;-6;32000000000;150000000000" "-9;-6;32000000000;150000000000" "-8;-6;32000000000;150000000000" "-7;-6;32000000000;150000000000" "-6;-6;32000000000;150000000000" "-5;-6;32000000000;150000000000" "-4;-6;32000000000;150000000000" "-3;-6;32000000000;150000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario116
cp -rf $TMPDIR/results/* $SLURM_SUBMIT_DIR/simulateScenario116
rm -rf $TMPDIR/*

#End of script

