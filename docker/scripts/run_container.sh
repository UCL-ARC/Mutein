#!/bin/bash

#
# run the container using the default "bash" command
#

set -eu

source settings

docker run -it \
    --volume ${HOST_DATA_FOLDER}:${GUEST_DATA_FOLDER} \
    mutein_tre
