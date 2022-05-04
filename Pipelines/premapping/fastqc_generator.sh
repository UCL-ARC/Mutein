#!/bin/bash

set -eu

source ~/.mutein_settings

for DATASET in $(cat datasets/active_datasets)
do
    for ACCESSION in $(cat datasets/${DATASET}/active_accessions)
    do
        echo "fastqc -o datasets/${DATASET}/${ACCESSION}_fastqc -t 4 datasets/${DATASET}/${ACCESSION}/*.fastq"
    done
done \
> fastqc_joblist

TOTAL_TASKS=$(cat fastqc_joblist | wc --lines)

echo qsub -t 1-${TOTAL_TASKS} ${MUT_DIR}/Pipelines/premapping/fastqc_runner.sh
