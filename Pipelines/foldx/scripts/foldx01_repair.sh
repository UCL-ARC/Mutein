#!/bin/bash -l

# Batch script to run an array job under.

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=2:00:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=3G

# Request 15 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=15G

# Set the name of the job.
#$ -N foldx-repair

# Email myself the job status
#$ -m be

# Set the working directory to somewhere in your scratch space.  
#  This is a necessary step as compute nodes cannot write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID.
## $ -wd /home/ucbtlcr/Scratch/workspace

# Load the necessary python libraries
module load python3/recommended
module load foldx

# get the script inputs
pdb=$1
jobname=$2

# cd /home/ucbtlcr/Project549/Project549_02/nodes
cd $PWD
python foldx01_repair.py $pdb $jobname R MYR
