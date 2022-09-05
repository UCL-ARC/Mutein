#!/bin/bash -l
#
# Batch script for the Mutein pipeline - array jobs
#
# inputs that are overridden from pipeline script
#$ -l h_rt=5:00:0
#$ -wd /home/xxx/Scratch/workspace
#$ -N genes-proteins
#
# inputs in script only
#$ -l mem=4G
#$ -l tmpfs=15G
#
# DONT Email myself the job status - too many, outlook told me off
#
# First we want to make sure we know what params have been passed to the script
echo "~~~~~~~~~~~~~~~~~~ Parameters passed to the script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
# get the script inputs
script=$0
inputs=$1
pyscript=$2
workspace=$3
fullcall=$4

echo "$fullcall"
echo "SCRIPT=$script"
echo "PYSCRIPT=$pyscript"
echo "INPUTS=$inputs"
echo "TASKNUM=$SGE_TASK_ID"

cd $workspace
echo "CURRENT=$PWD"
# Load the necessary python libraries
module load python3/recommended
module load foldx
echo "~~~~~~~~~~~~~~~~~~ Call python script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python "$pyscript" "$inputs@task=$SGE_TASK_ID"
echo "MUTEIN SCRIPT ENDED"

