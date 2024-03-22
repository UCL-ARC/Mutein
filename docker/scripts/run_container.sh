#!/bin/bash

#
# test run the container on the build machine and activate the default "bash" interactive command prompt
#

set -eu

source settings

docker run -it \
    --volume ${HOST_DATA_FOLDER}:${GUEST_DATA_FOLDER} \
    mutein_tre
