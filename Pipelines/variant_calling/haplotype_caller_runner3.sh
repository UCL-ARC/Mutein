#!/bin/bash -l

#$ -l h_rt=3:00:0
#$ -l mem=9G
#$ -l tmpfs=1G
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -o sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -N gengvcf
###$ -tc 10

set -eu
source ~/.mutein_settings
conda activate ${MUT_PREFIX}gatk4

JOBLIST=$1

head -n ${SGE_TASK_ID} ${JOBLIST} | tail -n 1 | bash
