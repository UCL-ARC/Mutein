#!/bin/bash

set -eu
source ~/.mutein_settings

for DATASET in $(cat datasets/active_datasets)
do
    for ACCESSION in $(cat datasets/${DATASET}/active_accessions)
    do
        REF=reference/${MUT_REFERENCE}
        BAM=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_aln_sort.bam
        VCF=datasets/${DATASET}/${ACCESSION}/${ACCESSION}_haplo.vcf

        #use -L to specify target regions, eg for exome/bait capture
        echo -n "gatk -T HaplotypeCaller -R ${REF} -I ${BAM} -o ${VCF} "
        echo    "--genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30"
    done
done \
> haplo_joblist

TOTAL_TASKS=$(cat haplo_joblist | wc --lines)

echo qsub -t 1-${TOTAL_TASKS} ${MUT_DIR}/Pipelines/variant_calling/haplotype_caller_runner.sh
