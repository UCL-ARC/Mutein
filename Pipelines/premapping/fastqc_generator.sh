#!/bin/bash

set -eu

source ~/.mutein_settings

if [ $# -eq 0 ]; then
    JOBLIST_BASE=haplo_joblist
else
    JOBLIST_BASE=haplo_joblist_$1
fi

rm -f ${JOBLIST_BASE}

for DATASET in $(cat datasets/active_datasets)
do
    for ACCESSION in $(cat datasets/${DATASET}/active_accessions)
    do
        echo "fastqc -t 4 datasets/${DATASET}/${ACCESSION}/*.fastq.gz" >> ${JOBLIST_BASE}
    done
done

TOTAL_TASKS=$(cat ${JOBLIST_BASE} | wc --lines)

#for now just print qsub command to console for user to manually paste and run
#but we could call it directly, eg wrapped it in a screen session,
#use -sync y to catch errors etc
#eg screen -S factqc -d -m qsub -sync y...
echo qsub -t 1-${TOTAL_TASKS} ${MUT_DIR}/Pipelines/premapping/fastqc_runner.sh
