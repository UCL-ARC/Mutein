#!/bin/bash -l

# Batch script to make the paramater file for the variants

# inputs that are overridden from pipeline script
#$ -l h_rt=0:10:0
#$ -wd /home/ucbtlcr/Scratch/workspace

# inputs that are in the script only
#$ -l mem=1G
#$ -l tmpfs=15G
#$ -N foldx-vparams
# Email myself the job status
#$ -m be

# Load the necessary python libraries
module load python3/recommended

pdb="pdb="$1
jobname="name="$2
variant="variant="$5
variantfile="variantfile="$6
repairs="repairs="$7

cd ~/MuteinPipeline/foldx/scripts/
python foldx05_vparams.py $pdb $jobname $variant $variantfile $repairs


