#!/bin/bash
#
# Batch script to turn the dataset genes into protein structures
#
# inputs that are overridden from pipeline script
#$ -l h_rt=5:00:0
#$ -wd /home/xxx/Scratch/workspace
#
# inputs in script only
#$ -l mem=3G
#$ -l tmpfs=15G
#$ -N genes-proteins
# Email myself the job status
#$ -m be
#
# First we want to make sure we know what params have been passed to the script
echo "~~~~~~~~~~~~~~~~~~ Parameters passed to the script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
# get the script inputs
script=$0
inputs=$1
pyscript="$2"
echo "SCRIPT=$script"
echo "WORKDIR=$work_dir"
echo "INPUTS=$inputs"
echo "PYSCRIPT=$pyscript"
echo "CURRENTDIR=$PWD"
#
# Load the necessary python libraries
module load python3/recommended
echo "~~~~~~~~~~~~~~~~~~ Call python script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python $pyscript $inputs