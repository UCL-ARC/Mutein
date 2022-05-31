#!/bin/bash

#$ -l h_rt=4:00:0
#$ -l mem=6G
#$ -l tmpfs=6G
###$ -pe smp 1
#$ -cwd
#$ -V
#$ -o sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e sge_logs/$JOB_NAME.$JOB_ID.$TASK_ID.err
###$ -tc 10

set -eu
source ~/.mutein_settings
conda activate ${MUT_PREFIX}bwa

#extract gzipped file
gunzip reference/${MUT_REFERENCE}

#recompress with bgzip
bgzip reference/${MUT_REFERENCE/\.gz/}

#run approx ??? hour
samtools faidx reference/${MUT_REFERENCE}

#run time approx 1 hour
bwa index -a bwtsw reference/${MUT_REFERENCE}
conda deactivate ${MUT_PREFIX}bwa

conda activate ${MUT_PREFIX}gatk4

rm -f reference/${MUT_REFERENCE}.dict
gatk CreateSequenceDictionary -R reference/${MUT_REFERENCE}

conda deactivate ${MUT_PREFIX}gatk4
