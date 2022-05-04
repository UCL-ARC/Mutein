#!/usr/bin/env bash

#
# downloads the Keogh et al 2018 dataset
# you must first manually download the metadata and accession list file as described below
# save these two files under datasets/keogh2018 before running this script
# currently requires: bash, enaBrowserTools
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
source ~/.mutein_settings

#may become command line arguments in future
DATA_DIR=datasets/keogh2018
ACC_FILE=SRR_Acc_List.txt

module load ${MUT_CONDA_MODULE}

cd "${DATA_DIR}"

for x in $(cat "${ACC_FILE}")
do
    enaDataGet -f fastq ${x}
done \
> download-${TIMESTAMP}.stdout \
2> download-${TIMESTAMP}.stderr

#set all samples as active (manually edit active_acc_list to work selectly)
cp ${ACC_FILE} active_acc_list
