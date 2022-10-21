#!/bin/bash -l

#
# this run inside the container upon start up
# see the separate def file for how to container was originally created
#

echo start setup_container
modprobe fuse
#mkdir -p /tmp/mutein
#cat config/encfs_password.tmp | encfs --standard --stdinpass $1 /tmp/mutein
#echo after encfs

#launch interactive session
bash
