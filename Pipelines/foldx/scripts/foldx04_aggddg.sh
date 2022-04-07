#!/bin/bash -l

# Batch script to run aggregate the background mutations

# inputs that are overridden from pipeline script
#$ -l h_rt=0:10:0
#$ -wd /home/ucbtlcr/Scratch/workspace

# inputs that are in the script only
#$ -l mem=1G
#$ -l tmpfs=15G
#$ -N foldx-aggregate
# Email myself the job status
#$ -m be

# Load the necessary python libraries
module load python3/recommended

pdb="pdb="$1
jobname="name="$2
repairs="repairs="$7

cd ~/MuteinPipeline/foldx/scripts/
python foldx04_aggddg.py $pdb $jobname $repairs


