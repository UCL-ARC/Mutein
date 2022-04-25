#!/bin/bash -l

# Batch script to turn the dataset genes into protein structures

# inputs that are overridden from pipeline script
#$ -l h_rt=5:00:0
#$ -wd /home/xxx/Scratch/workspace

# inputs in script only
#$ -l mem=3G
#$ -l tmpfs=15G
#$ -N genes-proteins
# Email myself the job status
#$ -m be

# Load the necessary python libraries
module load python3/recommended
# pip install bioservices

# get the script inputs
work_dir=$1
echo $work_dir
inputs=$2
echo $inputs
# cd ~/MuteinPipeline/geneprot/scripts/
cd $work_dir
python pipeline02_genestoproteins.py $inputs
