#!/bin/bash

#
# generate joblist file for germline-type variant calling with gatk4 haplotypecaller
#

set -eu
source ~/.mutein_settings

if [ $# -eq 0 ]; then
    JOBLIST_BASE=haplo_joblist
else
    JOBLIST_BASE=haplo_joblist_$1
fi

rm -f ${JOBLIST_BASE}

SGE_TASK_ID=0
for DATASET in $(cat datasets/active_datasets)
do
    for ACCESSION in $(cat datasets/${DATASET}/active_accessions)
    do
        SGE_TASK_ID=$((SGE_TASK_ID+1))
        REF=reference/${MUT_REFERENCE}
        BAM=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_aln_sort.bam
        VCF=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_haplo.vcf
        INT=datasets/${DATASET}/${ACCESSION}/intervals.list
        ln -s ../intervals.list ${INT} #make all samples use the same intervals list for now

        #https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
        #use -L to specify target regions, eg for exome/bait capture
        #https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists
        echo -n "gatk --java-options \"-Xmx8g\" HaplotypeCaller " >> ${JOBLIST_BASE}
        echo -n "--native-pair-hmm-threads 4 -L ${INT} -ip 10 "   >> ${JOBLIST_BASE}
        echo    "-R ${REF} -I ${BAM} -O ${VCF} -ERC GVCF"         >> ${JOBLIST_BASE}
    done
done

TOTAL_TASKS=$(cat ${JOBLIST_BASE} | wc --lines)

echo qsub ${MUT_DIR}/Pipelines/variant_calling/haplotype_caller_prep.sh
echo qsub -t 1-${TOTAL_TASKS} ${MUT_DIR}/Pipelines/variant_calling/haplotype_caller_runner.sh ${JOBLIST_BASE}
