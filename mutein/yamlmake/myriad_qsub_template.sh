#!/bin/bash -l
#$ -N {jobname}
#$ -cwd
#$ -V
#$ -t 1-{njobs}
{_tc} -tc {maxrun}
#$ -o {log_dir}/{jobname}.$TASK_ID.out
#$ -e {log_dir}/{jobname}.$TASK_ID.err
#$ -l h_rt={time}
#$ -l mem={mem}
#$ -l tmpfs={tmpfs}
{_pe} -pe {pe} {cores}
{python} {yamlmake} --qsub {jobfile}
