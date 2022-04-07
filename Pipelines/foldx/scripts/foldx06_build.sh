#!/bin/bash -l

# Batch script to run an array of variant processes for ddg

# inputs that are overridden from pipeline script
#$ -l h_rt=3:00:0
#$ -t 1-2
#$ -wd /home/ucbtlcr/Scratch/workspace

# inputs that are in the script only
#$ -l mem=4G
#$ -l tmpfs=15G
#$ -N foldx-build
# Email myself the job status
#$ -m be

# Load the necessary python libraries
module load python3/recommended
module load foldx

# Parse parameter file to get variables.
pdb=$1
jobname=$2
rows=$3
a='/home/'$USER'/MuteinPipeline/foldx/thruputs/'
b='/variant_params.txt'
d=$1
paramfile=${a}${pdb}${b}

number=$SGE_TASK_ID

ipdb="`sed -n ${number}p $paramfile | awk '{print $1}'`"
ichain="`sed -n ${number}p $paramfile | awk '{print $2}'`"
imutation="`sed -n ${number}p $paramfile | awk '{print $3}'`"
irow="`sed -n ${number}p $paramfile | awk '{print $4}'`"

jobnamex="name="$2
pdb="pdb="$ipdb
rw="row="$irow
mutation="mutation="$imutation
repairs="repairs="$7


cd ~/MuteinPipeline/foldx/scripts/
python foldx06_build.py $pdb $jobnamex $mutation $rw $repairs

