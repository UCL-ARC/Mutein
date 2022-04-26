#!/bin/bash -l

# Batch script to run a repair job on pdb files

# inputs that are overridden from pipeline script
#$ -l h_rt=5:00:0
#$ -wd /home/xxx/Scratch/workspace

# inputs in script only
#$ -l mem=3G
#$ -l tmpfs=15G
#$ -N foldx-repair
# Email myself the job status
#$ -m be

# First we want to make sure we know what params have been passed to the script
echo "~~~~~~~~~~~~~~~~~~ Parameters passed to the script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
for var in "$@"
do
    echo "$var"
done

# Load the necessary python libraries
module load python3/recommended
module load foldx

# get the script inputs
echo $0
work_dir=$1
echo "work directory = $work_dir"
inputs=$2
echo "inputs = $inputs"
cd $work_dir
python foldx01_repair.py $inputs
