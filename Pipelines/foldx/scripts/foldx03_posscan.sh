#!/bin/bash -l



# Batch script to run an array job under.

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=3:00:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=4G

# Request 15 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=15G

# Set up the job array.   #$ -t 1-$rows
#$ -t 1-2

# Set the name of the job.
#$ -N foldx-posscan

# Email myself the job status
#$ -m be

# Set the working directory to somewhere in your scratch space.  
#  This is a necessary step as compute nodes cannot write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID.
#$ -wd /home/ucbtlcr/Scratch/workspace

# Load the necessary python libraries
module load python3/recommended
module load foldx

# Parse parameter file to get variables.
jobname=$2
rows=$3
a='/home/ucbtlcr/MuteinPipeline/foldx/interim/'
b='/params.txt'
d=$1
paramfile=${a}${jobname}${b}

number=$SGE_TASK_ID

ipdb="`sed -n ${number}p $paramfile | awk '{print $1}'`"
ichain="`sed -n ${number}p $paramfile | awk '{print $2}'`"
imutation="`sed -n ${number}p $paramfile | awk '{print $3}'`"
arow="`sed -n ${number}p $paramfile | awk '{print $4}'`"

jobnamex="name="$2
pdb="pdb="$ipdb
mutation="mutation="$imutation
xrow = "row="$arow
newrowa = "rowa="${number}p
newrowb = "rowb="$SGE_TASK_ID


cd ~/MuteinPipeline/foldx/scripts/
python foldx03_posscan.py A=10 Z=15 $newrowa $newrowb $pdb $jobnamex $mutation $xrow

