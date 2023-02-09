#!/bin/bash -l
#$ -N {jobname}
#$ -cwd
#$ -V
#$ -t 1-{njobs}
#### -tc 20 #max concurrent if required
#$ -o {log_dir}/{jobname}.$TASK_ID.out
#$ -e {log_dir}/{jobname}.$TASK_ID.err
#$ -l h_rt={time}
#$ -l h_vmem={mem},tmem={mem}
#$ -l tscratch={tmpfs}
#$ -pe {pe} {cores}
{python} {yamlmake} --qsub {jobfile}
