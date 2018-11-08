#!/usr/bin/env bash
#SBATCH -p glenn
#SBATCH -A C3SE2018-1-15
#SBATCH -J simulateScenario217
#SBATCH -N 1
#SBATCH -t 0-02:00:00
#SBATCH -o simulateScenario217.stdout
#SBATCH -e simulateScenario217.stderr


module load matlab

cp  -r $SLURM_SUBMIT_DIR/* $TMPDIR
cd $TMPDIR

array=(  "2;7;64000000000;200000000000" "3;7;64000000000;200000000000" "4;7;64000000000;200000000000" "5;7;64000000000;200000000000" "6;7;64000000000;200000000000" "7;7;64000000000;200000000000" "8;7;64000000000;200000000000" "9;7;64000000000;200000000000" "10;7;64000000000;200000000000" "-10;8;64000000000;200000000000" "-9;8;64000000000;200000000000" "-8;8;64000000000;200000000000" "-7;8;64000000000;200000000000" "-6;8;64000000000;200000000000" "-5;8;64000000000;200000000000" "-4;8;64000000000;200000000000" )
for i in "${array[@]}"
do
    arr=(${i//;/ })
    echo ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}
    RunMatlab.sh -o "-nodesktop -nosplash -singleCompThread -r \"simulateScenario(${arr[0]},${arr[1]},${arr[2]},${arr[3]});\"" & 
    sleep 0.1
done

wait

mkdir $SLURM_SUBMIT_DIR/simulateScenario217
cp -rf $TMPDIR/* $SLURM_SUBMIT_DIR/simulateScenario217
rm -rf $TMPDIR/*

#End of script

