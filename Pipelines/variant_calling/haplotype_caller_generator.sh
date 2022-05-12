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

        #https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
        #use -L to specify target regions, eg for exome/bait capture
        echo -n "gatk --java-options \"-Xmx4g\" HaplotypeCaller "
        echo    "-R ${REF} -I ${BAM} -O ${VCF} -ERC GVCF"
    done
done \
> haplo_joblist

TOTAL_TASKS=$(cat haplo_joblist | wc --lines)

echo qsub ${MUT_DIR}/Pipelines/variant_calling/haplotype_caller_prep.sh
echo qsub -t 1-${TOTAL_TASKS} ${MUT_DIR}/Pipelines/variant_calling/haplotype_caller_runner.sh
