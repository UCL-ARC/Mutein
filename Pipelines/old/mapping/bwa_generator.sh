#!/bin/bash

#
# generate joblist files and qsub commands for bwa indexing, mapping, sorting  and bam indexing
#

set -eu
set -o pipefail
source ~/.mutein_settings

#optional basename for joblist files
if [ $# -eq 0 ]
then
    JOBLIST_BASE=bwa_joblist
else
    JOBLIST_BASE=$1
fi

JOBLIST1=${JOBLIST_BASE}_step1
JOBLIST2=${JOBLIST_BASE}_step2

rm -f ${JOBLIST1} ${JOBLIST2}

SGE_TASK_ID=0
for DATASET in $(cat datasets/active_datasets)
do
    for SUBSET in $(cat datasets/${DATASET}/active_subsets)
    do
        for ACCESSION in $(cat datasets/${DATASET}/${SUBSET}/active_accessions)
        do
            SGE_TASK_ID=$((SGE_TASK_ID+1))
            READS1=datasets/${DATASET}/${SUBSET}/${ACCESSION}/${ACCESSION}_1.fastq.gz
            READS2=datasets/${DATASET}/${SUBSET}/${ACCESSION}/${ACCESSION}_2.fastq.gz
            REF=reference/${MUT_REFERENCE}
            BAMOUT=datasets/${DATASET}/${SUBSET}/${ACCESSION}/${ACCESSION}_aln_sort.bam

            #array job commands for step1:
            #map with bwa mem, add readgroup tag to say "these are all from the same sample"
            #sort mapped read by mapping position, output to bam file
            echo -n "bwa mem -t 4 ${REF} ${READS1} ${READS2} | "              >> ${JOBLIST1}
            echo -n "samtools addreplacerg -r \"ID:${DATASET}_${SUBSET}_${ACCESSION}\"" >> ${JOBLIST1}
            echo -n " -r \"SM:${DATASET}_${SUBSET}_${ACCESSION}\" - | "                 >> ${JOBLIST1}
            echo    "samtools sort -T /tmp/${USER}-bwa-${SGE_TASK_ID} -O bam -m 2G - > ${BAMOUT}"  >> ${JOBLIST1}

            #array job commands for step2:
            #index the sorted bam file
            echo "samtools index ${BAMOUT}"              >> ${JOBLIST2}
        done
    done
done

JOBNAME=bwa-index-$(mutein_random_id)
echo qsub -N ${JOBNAME} ${MUT_DIR}/Pipelines/mapping/bwa_index.sh
echo capture_qacct.sh ${JOBNAME}

JOBNAME=bwa-step1-$(mutein_random_id)
TOTAL_TASKS1=$(cat ${JOBLIST1} | wc --lines)
echo qsub -N ${JOBNAME} -t 1-${TOTAL_TASKS1} ${MUT_DIR}/Pipelines/mapping/bwa_runner_step1.sh ${JOBLIST1}
echo capture_qacct.sh ${JOBNAME}

JOBNAME=bwa-step2-$(mutein_random_id)
TOTAL_TASKS2=$(cat ${JOBLIST2} | wc --lines)
echo qsub -N ${JOBNAME} -t 1-${TOTAL_TASKS2} ${MUT_DIR}/Pipelines/mapping/bwa_runner_step2.sh ${JOBLIST2}
echo capture_qacct.sh ${JOBNAME}
