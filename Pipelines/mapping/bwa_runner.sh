#!/bin/bash -l

#$ -l h_rt=3:00:0
#$ -l mem=2G
#$ -l tmpfs=4G
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -o sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -N bwa
###$ -tc 10

set -eu
source ~/.mutein_settings
module load ${MUT_CONDA_MODULE}
conda activate bwa

mkdir -p /tmp/${USER}-sam

head -n ${SGE_TASK_ID} bwa_joblist | tail -n 1 | bash
