#!/bin/bash

#
# download the reference genome sequence analysis set FASTA file from NCBI
# currently requires: bash, wget, md5sum
# currently downloads to a subfolder of the current working directory
# see the README_analysis_sets.txt file for explanation of what the analysis set includes
#

#quit on error or undefined variable
set -eu

SRC_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
REF_FILE=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
MD5_FILE=md5checksums.txt
DOC_FILE=README_analysis_sets.txt

DST_FOLDER=reference
CHK_FILE=md5_GRCh38_no_alt_plus_hs38d1.txt

#create folder to receive the reference sequence if not already present
mkdir -p ${DST_FOLDER} && cd ${DST_FOLDER}

#download the actual files
wget ${SRC_URL}/${MD5_FILE}
wget ${SRC_URL}/${DOC_FILE}
wget ${SRC_URL}/${REF_FILE}

#extract just the relevant checksums
cat ${MD5_FILE} | grep -e "${REF_FILE}" > ${CHK_FILE}

#perform the check
md5sum --check ${CHK_FILE}
