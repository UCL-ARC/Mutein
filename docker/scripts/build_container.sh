#!/bin/bash

#
# script used to build the docker container image on the docker build VM
# run from the directory containing the Dockerfile

set -eu

usage="usage: ./scripts/build_container.sh   #ie run from the 'docker' folder"

if [[ "$#" != "0" ]] ; then
    echo "$usage"
    exit
fi

#copy the github keys into the build context
mkdir -p secrets
chmod 0700 secrets
cp /vagrant/untracked/id_ed25519 secrets/id_ed25519
chmod 0600 secrets/id_ed25519
cp /vagrant/untracked/id_ed25519.pub secrets/id_ed25519.pub
chmod 0644 secrets/id_ed25519.pub

docker build -t mutein_tre .
docker image list