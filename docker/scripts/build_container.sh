#!/bin/bash

#
# script used to build the docker container image on the docker build VM
# run from the directory containing the Dockerfile
#

set -eu

usage="usage: ./scripts/build_container.sh   #ie run from the 'docker' folder"

if [[ "$#" != "0" ]] ; then
    echo "$usage"
    exit
fi

docker build -t mutein_tre .
docker image list