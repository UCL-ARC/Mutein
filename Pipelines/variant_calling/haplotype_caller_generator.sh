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

JOBLIST1=${JOBLIST_BASE}_step1
JOBLIST2=${JOBLIST_BASE}_step2
JOBLIST3=${JOBLIST_BASE}_step3

rm -f ${JOBLIST1} ${JOBLIST2} ${JOBLIST3}

SGE_TASK_ID=0
for DATASET in $(cat datasets/active_datasets)
do
    GDB=datasets/${DATASET}/gatk_db
    FINALVCF=datasets/${DATASET}/${DATASET}.vcf
    
    #https://gatk.broadinstitute.org/hc/en-us/articles/360035889971
    echo -n "rm -rf ${GDB} && gatk GenomicsDBImport" >> ${JOBLIST2}

    for ACCESSION in $(cat datasets/${DATASET}/active_accessions)
    do
        SGE_TASK_ID=$((SGE_TASK_ID+1))
        REF=reference/${MUT_REFERENCE}
        BAM=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_aln_sort.bam
        VCF=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_haplo.vcf
        INT=datasets/${DATASET}/${ACCESSION}/intervals.list

        rm -f ${INT}
        ln -s ../intervals.list ${INT} #make all samples use the same intervals list for now

        #https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
        #use -L to specify target regions, eg for exome/bait capture
        #https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists
        echo -n "gatk --java-options \"-Xmx8g\" HaplotypeCaller " >> ${JOBLIST1}
        echo -n "--native-pair-hmm-threads 4 -L ${INT} -ip 10 "   >> ${JOBLIST1}
        echo    "-R ${REF} -I ${BAM} -O ${VCF} -ERC GVCF"         >> ${JOBLIST1}

        echo -n " -V ${VCF}" >> ${JOBLIST2}
    done

    echo " --genomicsdb-workspace-path ${GDB} -L datasets/${DATASET}/intervals.list -ip 10" >> ${JOBLIST2}

    echo "gatk GenotypeGVCFs -R ${REF} -V gendb://${GDB} -O ${FINALVCF}" >> ${JOBLIST3}

done

TOTAL_TASKS1=$(cat ${JOBLIST1} | wc --lines)
TOTAL_TASKS2=$(cat ${JOBLIST2} | wc --lines)
TOTAL_TASKS3=$(cat ${JOBLIST3} | wc --lines)

echo qsub ${MUT_DIR}/Pipelines/variant_calling/haplotype_caller_prep.sh

echo qsub -t 1-${TOTAL_TASKS1} ${MUT_DIR}/Pipelines/variant_calling/haplotype_caller_runner1.sh ${JOBLIST1}
echo qsub -t 1-${TOTAL_TASKS2} ${MUT_DIR}/Pipelines/variant_calling/haplotype_caller_runner2.sh ${JOBLIST2}
echo qsub -t 1-${TOTAL_TASKS3} ${MUT_DIR}/Pipelines/variant_calling/haplotype_caller_runner3.sh ${JOBLIST3}
