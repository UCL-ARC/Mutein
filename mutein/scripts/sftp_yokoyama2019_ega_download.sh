#!/bin/bash

#
# download yokoyama from ega outbox using sftp
#
# -r recursively copy
# -a attempt to continue interrupted transfers
# -F path to ssh config file
#
# to be run in /SAN/medic/Mutein/549_mutein/datasets/yokoyama2019
#

timestamp=$(date +%Y%m%d-%H%M%S)

sftp -ra -F ./ssh_config ega-outbox 2> ${timestamp}.err > ${timestamp}.out <<EOF
get EGAD00001004462
get EGAD00001004464
bye
EOF


