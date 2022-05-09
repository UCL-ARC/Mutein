#!/bin/bash

set -eu

source ~/.mutein_settings

for DATASET in $(cat datasets/active_datasets)
do
    for ACCESSION in $(cat datasets/${DATASET}/active_accessions)
    do
        echo "fastqc -t 4 datasets/${DATASET}/${ACCESSION}/*.fastq.gz"
    done
done \
> fastqc_joblist

TOTAL_TASKS=$(cat fastqc_joblist | wc --lines)

#for now just print qsub command to console for user to manually paste and run
#but we could call it directly, eg wrapped it in a screen session,
#use -sync y to catch errors etc
#eg screen -S factqc -d -m qsub -sync y...
echo qsub -t 1-${TOTAL_TASKS} ${MUT_DIR}/Pipelines/premapping/fastqc_runner.sh
