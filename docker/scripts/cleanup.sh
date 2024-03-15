#!/bin/bash

set -eu

docker system prune -f --volumes
docker system prune -f --volumes
