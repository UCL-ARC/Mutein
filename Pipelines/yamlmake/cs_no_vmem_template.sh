#!/bin/bash -l
# see https://hpc.cs.ucl.ac.uk/wp-content/uploads/sites/21/2022/01/SGE-Guide_12_2021.pdf p19
# which notes that sometime h_vmem should be omitted as it causes memory errors
#$ -N {jobname}
#$ -cwd
#$ -V
#$ -t 1-{njobs}
#### -tc 10 #max concurrent if required
#$ -o {log_dir}/{jobname}.$TASK_ID.out
#$ -e {log_dir}/{jobname}.$TASK_ID.err
#$ -l h_rt={time}
#$ -l tmem={mem}
#$ -l tscratch={tmpfs}
#$ -pe {pe} {cores}
{python} {yamlmake} --qsub {jobfile}
