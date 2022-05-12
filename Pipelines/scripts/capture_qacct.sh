#!/bin/bash

#
# capture post job status info from qacct
# captures owner, taskid, "failed", "exit_status", "maxvmem" and "ru_wallclock"

JOB_NAME=$1
JOB_ID=$2

source ~/.mutein_settings

qacct -j "${JOB_ID}" \
| grep -e failed -e exit_status -e maxvmem -e ru_wallclock -e owner -e taskid -e jobnumber\
| paste - - - - - - - | grep -e ${USER} | grep -e "${JOB_ID}" \
> sge_logs/${JOB_NAME}.${JOB_ID}.qacct
