#!/bin/bash -l

# Batch script to create parameters for the background mutations

# inputs that are overridden from pipeline script
#$ -l h_rt=0:10:0
#$ -N foldx-params
#$ -wd /home/ucbtlcr/Scratch/workspace

# inputs in script only
#$ -l mem=1G
#$ -l tmpfs=15G
# Email myself the job status
#$ -m be

# Load the necessary python libraries
module load python3/recommended

pdb=$1
jobname=$2
chunk=$3

pdb="pdb="$1
jobname="name="$2
chunk="split="$3
repairs="repairs="$7

cd ~/MuteinPipeline/foldx/scripts/
python foldx02_makeparams.py $pdb $jobname $chunk $repairs
