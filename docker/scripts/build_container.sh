#!/bin/bash

#
# script used to build the docker container image on the docker build VM
# run from the directory containing the Dockerfile
#
# useful build options:
# --no-cache #to force full rebuild

set -eu

source settings

usage="usage: ./scripts/build_container.sh [--help|build-options]  #run from the 'docker' folder"

if [[ "$#" != "0" ]] ; then
    if [[ "$1" == "--help" || "$1" == "-h" ]] ; then
        echo "$usage"
        exit
    fi
fi

docker build $@ -t mutein_tre \
    --build-arg "username=${CONTAINER_USER}" \
    --build-arg "data_path=${GUEST_DATA_FOLDER}" \
    --build-arg "mamba_url=${MAMBA_URL}" \
    --build-arg "mamba_sh=${MAMBA_SH}" \
    .

#docker image list