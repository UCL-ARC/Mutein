#!/bin/bash

#
# capture post job status info from qacct
# jobid is not unique in the qacct database as the job counter seems to get reset
# and the database retains jobs from multiple counter "epochs"
# therefore we assume to jobname contains a UID and search by jobname instead
#

JOB_NAME=$1

source ~/.mutein_settings

qacct -j "${JOB_NAME}" \
| grep -e failed -e exit_status -e maxvmem -e ru_wallclock -e owner -e taskid -e jobnumber -e jobname\
| paste - - - - - - - - \
> sge_logs/${JOB_NAME}.qacct
