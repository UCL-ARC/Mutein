#!/usr/bin/env bash

#
# downloads the Yokoyama et al 2019 dataset
# currently requires: bash, pyega3
# run on headnode inside a screen session
# enter login details manually to prevent saving them to disk
#

set -eu
#source ~/.mutein_settings
export TIMESTAMP=$(date +%Y%m%d%H%M%S)

#may become command line arguments in future
DATA_DIR=datasets/yokoyama2019
ACC_FILE=EGA_acc_list #create manually

#module load ${MUT_CONDA_MODULE}
conda activate pyega3

mkdir -p "${DATA_DIR}"
cd "${DATA_DIR}"

for x in $(cat "${ACC_FILE}")
do
    pyega3 -ms 100000000000 fetch ${x}
done \
> download-${TIMESTAMP}.stdout \
2> download-${TIMESTAMP}.stderr

#set all samples as active (manually edit active_acc_list to work selectly)
#cp ${ACC_FILE} active_acc_list
