#!/bin/bash

#
# setup all the conda packages
#

set -eu

if [[ "$#" < "1" ]] ; then
    echo "sets up one or more conda environments from environment.yml type files"
    echo "example usage: setup_conda_envs.sh mutein/conda_envs/*.yml"
    exit
fi

for env_file in $@
do
    conda env create -f ${env_file} \
    || conda env update -f ${env_file} --prune
done
