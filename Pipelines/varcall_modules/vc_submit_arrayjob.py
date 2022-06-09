#!/usr/bin/env python

import sys
import os
import argparse
import subprocess

import varcall as vc

def count_tasks(filename,n_tasks):
    'count number of tasks in the file (number of lines excluding the header)'

    if filename == 'no-task-file':
        #n_tasks option must specify number of tasks if there is not taskfile
        assert n_tasks != None,'must provide either a taskfile or n_tasks'
        return int(n_tasks)

    f = open(filename)
    header = f.readline()
    assert header.startswith('['), "does not look like a task manifest file"
    count = 0
    for line in f: count += 1
    f.close()
    return count

def generate_environment(args,task_count):
    '''
    add parameters to the job's environment using the "preamble" part of the bash command
    '''

    env = os.environ.copy()
    env['vc_taskcount'] = str(task_count)
    for key,value in vars(args).items():
        if key in ['conf','jobname']:
            #don't pass these through at all
            continue
        elif key in vc.get_standard_args()+['taskfile','selectfile']:
            #prefix these resource request values with vc_
            #("time" etc seems a bit in risk of a name collision otherwise)
            env["vc_"+key] = value
        elif type(value) == str:
            #pass everything else through as is provided the value is a simple string
            env[key] = value
    return env

def submit_arrayjob():
    '''
    get command from last line of file, insert an up to date timestamp
    submit the job to qsub
    '''
    args = parse_args()

    task_count = count_tasks(args.taskfile,args.n_tasks)
    jobscript = os.path.join(os.environ['MUT_TEMPLATE_DIR'],args.jobname+'.qsub')
    assert os.path.exists(jobscript)

    env = generate_environment(args,task_count)

    cmd   = "qsub -cwd -V"

    if args.logs != "none":
        cmd += " -o '{args.logs}/$JOB_NAME.$TASK_ID.o' -e '{args.logs}/$JOB_NAME.$TASK_ID.e'"

    #add unique timestamp based id to jobname to avoid relying on JOBID
    jobname_uid = vc.unique_id()

    cmd += " -N {args.jobname}-{jobname_uid} -t 1-{task_count}"
    cmd += " -l h_rt={args.time} -l mem={args.mem} -l tmpfs={args.tmpfs} -pe smp {args.cores}"
    cmd += " {jobscript}"

    #update the timestamp
    cmd = cmd.format(**locals())

    print('\n'.join([key+'='+value for key,value in env.items()]) + '\n')
    print(cmd)

    #issue the command, ensure no immediate error encountered
    #subprocess.run(cmd.split(),check=True,env=env)

def parse_args():
    'get the filename of the array job specification file'
    parser = argparse.ArgumentParser(description="Submit a qsub array job defined by an arrayjob specification file previously made by a generator")
    parser.add_argument('--jobname',  required=True, type=str, help='Name of job and also jobscript template .qsub file')
    parser.add_argument('--n_tasks',  type=int, default=None, help='Specify number of tasks where no task file required')
    parser.add_argument('--taskfile',  default='no-task-file', type=str, help='Path to the task manifest file from generator function')
    parser.add_argument('--selectfile', default='no-selector-file', type=str, help='Path to the task select file to allow running only a subset of tasks')
    vc.add_conf_arg(parser)      #add --conf option
    vc.add_standard_args(parser) #add grid engine related arguments
    args = vc.parse_and_load_conf(parser) 

    return args

if __name__ == "__main__":
    globals()["submit_arrayjob"]()
