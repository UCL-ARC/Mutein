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
    --build-arg "REPO_FOLDER=${REPO_FOLDER}" \
    --build-arg "CONTAINER_USER=${CONTAINER_USER}" \
    --build-arg "GUEST_DATA_FOLDER=${GUEST_DATA_FOLDER}" \
    --build-arg "MAMBA_URL=${MAMBA_URL}" \
    --build-arg "MAMBA_SH=${MAMBA_SH}" \
    .

#docker image list