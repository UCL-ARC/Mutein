#!/bin/bash

#
# script used to build the docker container image on the docker build VM
# run from the directory containing the Dockerfile
#

set -eu

source settings

usage="usage: ./scripts/build_container.sh [--no-cache]  #run from the 'docker' folder"

nocache=""

if [[ "$#" != "0" ]] ; then
    if [[ "$1" == "--no-cache" ]] ; then
        nocache="$1"
    else
        echo "$usage"
        exit
    fi
fi

docker build "${nocache}" -t mutein_tre \
    --build-arg username="${CONTAINER_USER}" \
    --build-arg data_path="${GUEST_DATA_FOLDER}" \
    --build-arg mamba_url="${MAMBA_URL}" \
    --build-arg mamba_sh="${MAMBA_SH}" \
    .

docker image list