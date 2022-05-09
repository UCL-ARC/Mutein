#!/bin/bash -l

#$ -l h_rt=3:00:0
#$ -l mem=6G
#$ -l tmpfs=6G
###$ -pe smp 1
#$ -cwd
#$ -V
#$ -o sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -N bwa-index
###$ -tc 10

set -eu
source ~/.mutein_settings
module load ${MUT_CONDA_MODULE}
conda activate bwa

#run time approx 1 hour
bwa index -a bwtsw reference/${MUT_REFERENCE}
