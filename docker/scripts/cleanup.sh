#!/bin/bash

#
# cleanup cache docker volumes and images etc
#

set -eu

docker system prune -f --volumes
docker system prune -f --volumes
