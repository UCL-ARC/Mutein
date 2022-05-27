#!/usr/bin/env bash

#
# downloads the Yokoyama et al 2019 dataset
# currently requires: bash, pyega3
# run on headnode inside a screen session
# enter login details manually to prevent saving them to disk
#

set -eu
source ~/.mutein_settings
conda activate ${MUT_PREFIX}pyega3

#may become command line arguments in future
DATA_DIR=datasets/yokoyama2019
ACC_FILE=accession_list

#create accession list file, each item is a data subset contains multiple files
cat <<EOF > ${ACC_FILE}
EGAD00001004464
EGAD00001004462
EGAD00010001631
EOF

mkdir -p "${DATA_DIR}"
cd "${DATA_DIR}"

for x in $(cat "${ACC_FILE}")
do
    mkdir -p ${x} && cd ${x}
    pyega3 -ms 100000000000 fetch ${x}
    cd ..
done \
> download-${TIMESTAMP}.stdout \
2> download-${TIMESTAMP}.stderr

#set all samples as active (manually edit active_acc_list to work selectly)
cp ${ACC_FILE} active_acc_list
