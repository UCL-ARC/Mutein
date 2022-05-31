#!/bin/bash

#$ -l h_rt=0:30:0
#$ -l mem=1G
#$ -l tmpfs=15G
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -o sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.err
###$ -tc 10

set -eu
source ~/.mutein_settings
conda activate ${MUT_PREFIX}trim-galore
JOBLIST=$1

head -n ${SGE_TASK_ID} ${JOBLIST} | tail -n 1 | bash
