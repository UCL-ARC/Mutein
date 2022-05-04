#!/bin/bash -l

#$ -l h_rt=3:00:0
#$ -l mem=1G
#$ -l tmpfs=15G
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -o sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -N fastqc
#$ -tc 10

set -eu
source ~/.mutein_settings
module load ${MUT_CONDA_MODULE}
conda activate trim-galore

head -n ${SGE_TASK_ID} fastqc_joblist | tail -n 1 | bash
