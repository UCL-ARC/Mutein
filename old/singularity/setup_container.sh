#!/bin/bash -l

#
# this runs inside the container upon start up
# see the separate def file for how the container was originally created
# this script could probably be better implemented as a %startscript in the def file
#

echo start setup_container
modprobe fuse
#mkdir -p /tmp/mutein
#cat config/encfs_password.tmp | encfs --standard --stdinpass $1 /tmp/mutein
#echo after encfs

#launch interactive session
bash
