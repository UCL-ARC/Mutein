#!/bin/bash

#
# downloads the Keogh et al 2018 dataset
# edirect can be used to download the list of run ids
# but the full metadata must currently be downloaded manually
# you must first manually download the metadata and accession list file as described below
# save these two files under datasets/keogh2018 before running this script
# currently requires: bash, enaBrowserTools, entrez-direct
#

# in a browser visit: https://www.ncbi.nlm.nih.gov/sra/?term=SRP159015
# click on "send results to run selector"
# click on Select=>Total=>Metadata and save the SraRunTable.txt file
# click on Select=>Total=>Accession List and save the SRR_Acc_List.txt file
# see also https://www.ncbi.nlm.nih.gov/books/NBK179288/

set -eu
source ~/.mutein_settings
conda activate ${MUT_PREFIX}enabrowsertools

DATA_DIR=datasets/keogh2018/SRP159015 #only one subset for this paper
ACC_FILE=accession_list

mkdir -p ${DATA_DIR}
cd "${DATA_DIR}"

#grab list of sra run accessions from NCBI (across the pond)
esearch -db sra -query SRP159015 \
| efetch -format xml \
| xtract -pattern RUN -element PRIMARY_ID \
| awk '{print $1}' \
> ${ACC_FILE}

#download each accession as a fastq file from ENA (this side of the pond)
for x in $(cat "${ACC_FILE}")
do
    enaDataGet -f fastq ${x}
done \
> download-${TIMESTAMP}.stdout \
2> download-${TIMESTAMP}.stderr

#set all samples as active (manually edit active_acc_list to work selectly)
cp ${ACC_FILE} active_acc_list
