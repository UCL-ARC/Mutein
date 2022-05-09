#!/bin/bash -l

#$ -l h_rt=3:00:0
#$ -l mem=4G
#$ -l tmpfs=1G
###$ -pe smp 1
#$ -cwd
#$ -V
#$ -o sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -N haplcallr
###$ -tc 10

set -eu
source ~/.mutein_settings
module load ${MUT_CONDA_MODULE}
set +eu
conda activate gatk
set -eu

head -n ${SGE_TASK_ID} haplo_joblist | tail -n 1 | bash
