#!/usr/bin/env bash

set -eu
source ~/.mutein_settings

SRC_SAMPLE=SRR7762481
SRC_BASE=datasets/keogh2018/SRP159015/${SRC_SAMPLE}/${SRC_SAMPLE}
SRC_LINES=4000
SRC_INTERVALS=datasets/keogh2018/SRP159015/intervals.list
DST_SAMPLE=sample1

mkdir -p datasets/test/subset1/${DST_SAMPLE}
echo subset1 > datasets/test/active_subsets
echo ${DST_SAMPLE} > datasets/test/subset1/active_accessions
cp ${SRC_INTERVALS} datasets/test/subset1/intervals.list

zcat ${SRC_BASE}_1.fastq.gz \
| head -n ${SRC_LINES} \
| gzip \
> datasets/test/subset1/${DST_SAMPLE}/${DST_SAMPLE}_1.fastq.gz

zcat ${SRC_BASE}_2.fastq.gz \
| head -n ${SRC_LINES} \
| gzip \
> datasets/test/subset1/${DST_SAMPLE}/${DST_SAMPLE}_2.fastq.gz
