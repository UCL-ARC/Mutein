#!/bin/bash

#
# generate joblist files and qsub commands for bwa indexing, mapping, sorting  and bam indexing
#

set -eu
set -o pipefail
source ~/.mutein_settings

#optional basename for joblist files
if [ $# -eq 0 ]; then
    JOBLIST_BASE=bwa_joblist
else
    JOBLIST_BASE=bwa_joblist_$1
fi

rm -f ${JOBLIST_BASE}_step1 ${JOBLIST_BASE}_step2

SGE_TASK_ID=0
for DATASET in $(cat datasets/active_datasets)
do
    for ACCESSION in $(cat datasets/${DATASET}/active_accessions)
    do
        SGE_TASK_ID=$((SGE_TASK_ID+1))
        READS1=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_1.fastq.gz
        READS2=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_2.fastq.gz
        REF=reference/${MUT_REFERENCE}
        BAMOUT=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_aln_sort.bam

        #array job commands for step1:
        #map with bwa mem, add readgroup tag to say "these are all from the same sample"
        #sort mapped read by mapping position, output to bam file
        echo -n "bwa mem -t 4 ${REF} ${READS1} ${READS2} | "              >> ${JOBLIST_BASE}_step1
        echo -n "samtools addreplacerg -r \"ID:${DATASET}_${ACCESSION}\"" >> ${JOBLIST_BASE}_step1
        echo -n " -r \"SM:${DATASET}_${ACCESSION}\" - | "                 >> ${JOBLIST_BASE}_step1
        echo    "samtools sort -T /tmp/${USER}-bwa-${SGE_TASK_ID} -O bam -m 2G - > ${BAMOUT}"  >> ${JOBLIST_BASE}_step1

        #array job commands for step2:
        #index the sorted bam file
        echo "samtools index ${BAMOUT}"              >> ${JOBLIST_BASE}_step2
    done
done

#print the required qsub commands to screen for manual copy-pasting
echo qsub ${MUT_DIR}/Pipelines/mapping/bwa_index.sh

TOTAL_TASKS1=$(cat ${JOBLIST_BASE}_step1 | wc --lines)
echo qsub -t 1-${TOTAL_TASKS1} ${MUT_DIR}/Pipelines/mapping/bwa_runner_step1.sh ${JOBLIST_BASE}_step1

TOTAL_TASKS2=$(cat ${JOBLIST_BASE}_step2 | wc --lines)
echo qsub -t 1-${TOTAL_TASKS2} ${MUT_DIR}/Pipelines/mapping/bwa_runner_step2.sh ${JOBLIST_BASE}_step2
