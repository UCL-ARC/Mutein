#!/bin/bash -l

# Batch script to run an array job under.

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=0:10:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=1G

# Request 15 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=15G

# Set the name of the job.
#$ -N foldx-params

# Email myself the job status
#$ -m be

# Set the working directory to somewhere in your scratch space.  
#  This is a necessary step as compute nodes cannot write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID.
#$ -wd /home/ucbtlcr/Scratch/workspace

# Load the necessary python libraries
module load python3/recommended

pdb=$1
jobname=$2
chunk=$3

cd ~/MuteinPipeline/foldx/scripts/
python foldx02_makeparams.py $pdb $jobname A $chunk
