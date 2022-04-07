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

# Load the necessary python libraries
module load python3/recommended
module load foldx

# get the script inputs
pdb="pdb="$1
jobname="name="$2
echo $pdb, $jobname
repairs="repairs="$7

cd ~/MuteinPipeline/foldx/scripts/
python foldx01_repair.py $pdb $jobname $repairs
