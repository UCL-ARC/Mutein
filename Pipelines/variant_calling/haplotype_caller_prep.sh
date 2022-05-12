#!/bin/bash -l

#$ -l h_rt=1:00:0
#$ -l mem=4G
#$ -l tmpfs=1G
###$ -pe smp 1
#$ -cwd
#$ -V
#$ -o sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -N haplprep
###$ -tc 10

set -eu
source ~/.mutein_settings
module load ${MUT_CONDA_MODULE}
conda activate gatk4

#this should really be done earlier, when the reference is downloaded
rm -f reference/${MUT_REFERENCE}.dict
gatk CreateSequenceDictionary -R reference/${MUT_REFERENCE}
