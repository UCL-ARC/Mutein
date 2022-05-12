#!/bin/bash

set -eu

source ~/.mutein_settings

for DATASET in $(cat datasets/active_datasets)
do
    for ACCESSION in $(cat datasets/${DATASET}/active_accessions)
    do
        READS1=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_1.fastq.gz
        READS2=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_2.fastq.gz
        REF=reference/${MUT_REFERENCE}
        BAMOUT=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_aln_sort.bam
        echo -n "bwa mem -t 4 ${REF} ${READS1} ${READS2} | "
        echo -n "samtools addreplacerg -r \"ID:${DATASET}_${ACCESSION}\" - | "
        echo    "samtools sort -T /tmp/${USER}-sam -O bam -m 2G - > ${BAMOUT}"
    done
done \
> bwa_joblist

TOTAL_TASKS=$(cat bwa_joblist | wc --lines)

echo qsub ${MUT_DIR}/Pipelines/mapping/bwa_index.sh
echo qsub -t 1-${TOTAL_TASKS} ${MUT_DIR}/Pipelines/mapping/bwa_runner.sh
