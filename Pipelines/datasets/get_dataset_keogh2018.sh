#!/bin/bash

#
# downloads the Keogh et al 2018 dataset
# you must first manually download the metadata and accession list file as described below
# save these two files under datasets/keogh2018 before running this script
# currently requires: bash, sratools
#

# in a browser visit: https://www.ncbi.nlm.nih.gov/sra/?term=SRP159015
# click on "send results to run selector"
# click on Select=>Total=>Metadata and save the SraRunTable.txt file
# click on Select=>Total=>Accession List and save the SRR_Acc_List.txt file

# in the metadata SraRunTable.txt file seems like:
# Compare to Supplementary Table 4
# Isolate == LibraryName == Sample Name:
# Snnn_xxx_ttt
# nnn = sample number
# xxx = 001,002,003 = Control,Alzheimers,Lewy Body (see disease column)
# ttt: A,B,C,D,E,BLD = Cerebellum, Entorhinal cortex, Frontal Cortex, Medulla, Cingulate, Blood (see Cell_type column)

set -eu

SRATOOLS_DIR=~/software/sratoolkit.3.0.0-centos_linux64/bin
DATA_DIR=datasets/keogh2018
ACC_FILE=SRR_Acc_List.txt

export PATH=${SRATOOLS_DIR}:${PATH}

cd "${DATA_DIR}"

#download compressed sra files using prefetch command
#for x in $(cat "${ACC_FILE}")
#do
#    #prefetch --progress --resume yes --verify yes --max-size u --output-directory "${x}" "${x}"
#    break
#done

#problem: this seems to strip off the original Illumina readnames containing machine,lane,flow cell
#https://github.com/ncbi/sra-tools/issues/130
#for x in $(cat "${ACC_FILE}")
#do
#    fastq-dump --split-3 --skip-technical --gzip --outdir . "${x}"
#    break
#done

#might be original to EBI rather than NCBI?
#https://www.ebi.ac.uk/ena/browser/text-search?query=SRP159015
for x in $(cat "${ACC_FILE}")
do
    wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR776/005/${x}/${x}_1.fastq.gz
    wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR776/005/${x}/${x}_2.fastq.gz
    break
done
