import yaml
import re
import copy
import json
import os
import glob
import datetime
import time
import shutil
import subprocess
import sys
import random

#regex patterns to match non nested {placeholders} and {$environment variables}
ph_regx = r'\{[^\$=].*?\}'
#ph_regx = r'\{.+?\}'
en_regx = r'\{\$.+?\}'
gl_regx = r'\{=.+?\}'

#keys not allowed in config
nonconfig_keys = ['input','output','name']

#keys not allowed or required in rules
nonrule_keys = []
rule_keys = ['name','input','output','shell']

#placeholder first characters with special meaning
#$ environment variable
#= file glob: creates separate jobs from the same rule
###* file glob: creates file list within single job: not yet implemented
reserved_chrs = ['$','=']       #,'*']

default_global_config =\
{
    'log_dir':'fakemake_logs',       #name of subfolder for logging
    'log_prefix':'fm',               #log file/job name  prefix: qsub doesn't like numerical names
    'remote_delay_secs':'10',        #wait this long after remote jobs incase of latency
    'exec':'local',                  #default execution environment
    'conda_setup_command':  '',      #bash command to setup conda
    'conda':'fm_default_env',        #conda environment to activate if none specified by the rule
    'stale_output_file':'ignore',    #ignore,delete,recycle (also applies to symlinks)
    'stale_output_dir':'ignore',     #ignore,delete,recycle
    'failed_output_file':'stale',    #delete,recycle,stale,ignore (also applies to symlinks)
    'failed_output_dir':'stale',     #delete,recycle,stale,ignore    
    'missing_parent_dir':'create',   #ignore,create
    'recycle_bin':'recycle_bin',     #name of recycle bin folder
    'job_count':'FM_NJOBS',          #env variable: how many jobs spawned by current rule
    'job_number':'FM_JOB_NUMBER',    #env variable: 1 based job numbering within the current rule
    'bash_prefix':'source ~/.bashrc\nset -euo pipefail\nset +o history',
    'time':'02:00:00',          #$ -l h_rt={args.time}
    'mem':'4G',                 #$ -l mem={args.mem}
    'tmpfs':'10G',              #$ -l tmpfs={args.tmpfs}
    'pe':'smp',                 #$ -pe smp {threads}
    'cores':'1',
    'qsub_template':'default'   #default or path to your own
    #'normalize_paths':'true',                #true,false
}

def sub_vars(config,extra=None):
    'substitute all simple placeholders or fail trying'

    if extra == None: extra = config

    counter = 10
    while True:
        counter -= 1
        changed = sub_pholders(config,extra,'',ph_regx)
        if not changed: break
        assert counter > 0, 'unable to resolve all placeholders'

def sub_environ(config):
    'substitute in any placeholders for environment variables'
    return sub_pholders(config,os.environ,'$',en_regx)

def sub_globs(config,extra):
    'substitute in any globbing placeholders from extra'
    return sub_pholders(config,extra,'=',gl_regx)

def sub_pholders(config,extra,prefix,regx):
    '''
    find all placeholders in config values that match regx
    substitute in the corresponding value from extra
    raises exception if not found
    '''

    changed = False
    for key,value in config.items():
        while True:
            m = re.search(regx,value)
            if m == None: break

            name = m.group(0)[1:-1]
            if len(prefix) > 0: assert name.startswith(prefix)
            value = value[:m.start(0)] + extra[name[len(prefix):]] + value[m.end(0):]
            changed = True
        config[key] = value
    return changed

def show(item,label,indent=2):
    print(label)
    print(json.dumps(item,sort_keys=False,indent=indent))
    print()

def check_cwd(config):
    if 'expected_working_dir' in config:
        if os.path.realpath(os.getcwd()) != config['expected_working_dir']:
            print("warning: current path is not the expected working directory")
            print(f"expecting: {config['expected_working_dir']}")
            print(f"but found {os.path.realpath(os.getcwd())}")

def process(pipeline,path,config=None):
    #initially set config to default values if none provided
    if config == None:
        config = {}
        update_config(config,copy.deepcopy(default_global_config))

    counter = 0

    while counter < len(pipeline):
        item = pipeline[counter]
        counter += 1
        assert type(item) == dict
        assert len(item) == 1
        item_type = list(item.keys())[0]

        if item_type == 'rule':
            #show(item[item_type],item_type)
            process_rule(config,item[item_type])

        elif item_type == 'config':
            #add new config to the existing one, overriding any shared keys
            update_config(config,item[item_type])

        elif item_type == 'include':
            #load and insert the yaml items in place of the include item
            new_pipeline,new_path = load_pipeline(item[item_type],path) #new_path ignored
            counter -= 1
            del pipeline[counter]
            pipeline = pipeline[:counter] + new_pipeline + pipeline[counter:]

        elif item_type == 'module':
            #process a nested pipeline without affecting the config of any
            #following items
            new_pipeline,new_path = load_pipeline(item[item_type],path)
            new_config = copy.deepcopy(config)
            process(new_pipeline,new_path,config=new_config)

        #show(config,"config")

def load_pipeline(path,parent_file):
    'path is relative to the parent path'

    parent_path = os.path.dirname(parent_file)
    full_path = os.path.join(parent_path,path)

    assert full_path != parent_file
    return parse_yaml(full_path),full_path

def find_duplicates(dict_list):
    all_keys = set()

    for config in dict_list:
        for key in config.keys():
            assert not key in all_keys, f"duplicate key {key}"
            all_keys.add(key)

def setup_rule(config,rule):
    assert type(rule) == dict

    #check all keys and values are simple strings
    #check for forbidden keys
    for key,value in rule.items():
        assert key not in nonrule_keys
        assert type(key) == str
        if key not in ['input','output']:
            #general rule options must be simple strings
            assert type(value) == str
        else:
            #input and output can be simple strings or dictionaries of strings
            if type(value) != str:
                assert type(value) == dict
                for k,v in value.items():
                    assert type(k) == str
                    assert type(v) == str
    
    #check for required keys
    for key in rule_keys:
        assert key in rule

    #separate out and canonicalise input, output and shell
    input,output,shell = split_rule(rule)

    #merge temporary copy of config into the rule
    #such that rule keys overwrite any matching config keys
    _config = copy.deepcopy(config)
    _config.update(rule)
    rule.update(_config)

    #verify that rule, input and output do not share any keys
    find_duplicates([rule,input,output])

    #substitute any environment variables
    sub_environ(rule)
    sub_environ(input)
    sub_environ(output)
    sub_environ(shell)

    #substitute any fakemake variables except for the shell dict
    #which contains input/output variables yet to be determined
    sub_vars(rule)
    sub_vars(input,extra=rule)
    sub_vars(output,extra=rule)

    return input,output,shell

def split_rule(rule):
    #separate input and output from rule
    input = rule['input']
    output = rule['output']
    shell = rule['shell']
    del rule['input']
    del rule['output']
    del rule['shell']

    #convert simple form input/output into dictionary form
    if type(input) == str: input = {'input':input}
    if type(output) == str: output = {'output':output}
    shell = {'shell':shell}

    return input,output,shell

def update_config(config,new_config):
    #check all keys and values are simple strings
    #check for forbidden keys
    for key,value in new_config.items():
        assert type(key) == str
        assert type(value) == str
        assert key not in nonconfig_keys

    #substitute any environment variables
    sub_environ(new_config)

    #substitute any fakemake variables
    sub_vars(new_config)

    #create log directory if not already present
    if 'log_dir' in new_config:
        if not os.path.exists(new_config['log_dir']):
            os.makedirs(new_config['log_dir'])

    config.update(new_config)

def parse_yaml(fname):
    'parse yaml into list of items'
    with open(fname) as f:
        result = yaml.safe_load(f)

    assert type(result) == list
    
    return result

def esc_regx(value):
    'escape a literal that is being inserted into a regex'

    value = value.replace('.','\.')

    return value

def generate_job_list(rule,input,output):
    '''
    glob the first input pattern to file names
    fill out remaining input and output patterns with the matched names
    '''

    primary_key = list(input.keys())[0]
    primary = input[primary_key]

    input_glob = ''
    input_regx = '^'
    ditems = {}
    prev_end = 0

    #convert fakemake placeholders into glob and regex query formats
    for m in re.finditer(gl_regx,primary):
        name = m.group(0)[1:-1]
        start = m.start(0)

        input_glob += primary[prev_end:start] + '*'
        input_regx += esc_regx(primary[prev_end:start])

        #{=placeholder} defines a set of separate jobs
        name = name[1:]
        if name not in ditems:
            ditems[name] = True
            input_regx += '(?P<' + name + '>[^/]+)'
        else:
            #back reference to previous job or list placeholder
            input_regx += '(?P=' + name + ')'

        prev_end = m.end(0)

    input_glob += primary[prev_end:]
    input_regx += esc_regx(primary[prev_end:]) + '$'

    job_list = []

    #find paths using iglob, match to placeholders using regex
    for path in glob.iglob(input_glob):
        m = re.fullmatch(input_regx,path)
        assert m is not None

        job_input = copy.deepcopy(input)
        job_output = copy.deepcopy(output)

        sub_globs(job_input,m.groupdict())
        sub_globs(job_output,m.groupdict())

        job_list.append({"input":job_input,"output":job_output})

    return job_list

def check_input_mtimes(input):
    #verify all input paths present
    all_inputs_present = True
    newest_mtime = 0
    for item,path in input.items():
        if not os.path.exists(path) or is_stale(path):
            all_inputs_present = False
            continue

        mtime = os.path.getmtime(path)
        if mtime > newest_mtime: newest_mtime = mtime

    #return newest mtime if all inputs present otherwise None
    if all_inputs_present == False: return None
    return newest_mtime

def check_output_mtimes(output):
    all_outputs_present = True
    oldest_mtime = time.time() + 31e6 #dummy time 1 year in the future
    for item,path in output.items():
        if not os.path.exists(path):
            all_outputs_present = False
            continue

        mtime = os.path.getmtime(path)
        if mtime < oldest_mtime: oldest_mtime = mtime

    #return oldest mtime if all outputs present otherwise None
    if all_outputs_present == False: return None
    return oldest_mtime

def generate_shell_commands(rule,job_list,shell):
    name = rule['name']
    for job_numb,job in enumerate(job_list):
        #check all inputs present, return newest mtime
        newest_input = check_input_mtimes(job['input'])

        #one or more inputs missing, job not runnable
        if newest_input == None:
            #flag job for removal from the list
            job_list[job_numb] = None
            continue

        #all inputs present
        #check if any outputs missing or older than newest input
        oldest_output = check_output_mtimes(job['output'])

        #all outputs present and newer than newest input: no need to run
        if oldest_output != None and oldest_output > newest_input:
            #flag job for removal from the list
            job_list[job_numb] = None
            continue

        #remove and stale output files/symlinks/directories
        handle_stale_outputs(rule,job['output'])

        #create any missing output *parent* directories
        if rule['missing_parent_dir'] == 'create':
            create_output_dirs(rule,job['output'])

    #filter out deleted jobs
    job_list = [job for job in job_list if job is not None]
    shell_list = []

    for job_numb,job in enumerate(job_list):
        #merge input and output filenames into rule variables
        config = copy.deepcopy(rule)
        config.update(job['input'])
        config.update(job['output'])

        #substitute remaining placeholders in shell command
        shell_final = copy.deepcopy(shell)
        sub_vars(shell_final,config)
        #show(shell_final,"shell")
        shell_list.append(shell_final["shell"])

    return job_list,shell_list

def is_stale(path):
    'return true if mtime is set to the special stale flag'

    if os.path.getmtime(path) < datetime.datetime.fromisoformat('1981-01-01').timestamp():
        return True
    else:
        return False

def make_stale(path):
    '''
    set mtime of the item to special value to mark it as stale
    without deleting it (incase the contents need to be inspected first)
    '''

    dt_old = datetime.datetime.fromisoformat('1980-01-01').timestamp()
    os.utime(path, (dt_old, dt_old))

def recycle_item(config,path):
    'move the path into the recycle bin'
    
    recycle_bin = config['recycle_bin']

    assert path != recycle_bin

    new_name = os.path.join(recycle_bin,os.path.basename(path))

    #make sure not to overwrite existing recycle bin file
    if os.path.exists(new_name):
        new_name += ".%d"%random.randrange(1000000)
        if os.path.exists(new_name):
            raise Exception(f"File {new_name} already in recycle bin")

    #move to recycle bin
    shutil.move(path,new_name)

def create_output_dirs(config,outputs):
    '''
    if the parent folder of any output is missing create it
    '''

    for item in outputs:
        path = outputs[item]
        parent = os.path.dirname(path)
        if os.path.exists(parent): continue
        os.makedirs(parent)

def handle_stale_outputs(config,outputs):
    '''
    the job is going to be rerun therefore any output already present
    is treated as stale and recycled or deleted
    '''

    for item in outputs:
        path = outputs[item]
        if not os.path.exists(path): continue

        if os.path.islink(path) or os.path.isfile(path):
            if config['stale_output_file'] == 'delete':
                os.remove(path)
            elif config['stale_output_file'] == 'recycle':
                recycle_item(config,path)
            elif config['stale_output_file'] == 'ignore':
                continue
            else:
                raise Exception(f'unknown option for stale_output_file: {config["stale_output_file"]}')

        elif os.path.isdir(path):
            if config['stale_output_dir'] == 'delete':
                shutil.rmtree(path)
            elif config['stale_output_dir'] == 'recycle':
                recycle_item(config,path)
            elif config['stale_output_dir'] == 'ignore':
                continue
            else:
                raise Exception(f'unknown option for stale_output_dir: {config["stale_output_dir"]}')
                
        else:
            raise Exception(f'unsupported output type {path}')

def handle_failed_outputs(config,outputs):
    '''
    the job has failed therefore any outputs present are untrustworthy
    and should be removed, recycled, flagged as old etc
    '''

    for item in outputs:
        path = outputs[item]
        #print("handle_failed_outputs",path)
        if not os.path.exists(path): continue

        if os.path.islink(path) or os.path.isfile(path):
            #print('file or symlink')
            if config['failed_output_file'] == 'delete':
                #print('delete')
                os.remove(path)
            elif config['failed_output_file'] == 'recycle':
                #print('recycle')
                recycle_item(config,path)
            elif config['failed_output_file'] == 'stale':
                #print('stale')
                make_stale(path)
            elif config['failed_output_file'] == 'ignore':
                #print('ignore')
                continue
            else:
                raise Exception(f'unknown option for failed_output_file: {config["failed_output_file"]}')

        elif os.path.isdir(path):
            if config['failed_output_dir'] == 'delete':
                shutil.rmtree(path)
            elif config['failed_output_dir'] == 'recycle':
                recycle_item(config,path)
            elif config['failed_output_dir'] == 'stale':
                make_stale(path)
            elif config['failed_output_dir'] == 'ignore':
                continue
            else:
                raise Exception(f'unknown option for failed_output_dir: {config["failed_output_dir"]}')
                
        else:
            raise Exception(f'unsupported output type {path}')

def generate_full_command(config,shell):
    'generate the bash commands to run before the job commands'

    cmd_list = \
    [
        config['bash_prefix'],
        config['conda_setup_command'],
        f'conda activate {config["conda"]}',
        shell,
    ]

    return '\n'.join(cmd_list)

def generate_job_environment(config,job_numb,njobs):
    'copy local environment with a few adjustments'

    env = copy.deepcopy(os.environ)
    env[ config['job_count'] ] = f'{njobs}'
    env[ config['job_number'] ] = f'{job_numb+1}' # convert from 0 to 1 based to match SGE_TASK_ID

    return env

def write_jobfile(rule,shell_list):
    'save rule config and shell_list as json'

    fnamebase = f'{rule["log_prefix"]}{timestamp_now()}.{rule["name"]}'
    jobfile = os.path.join(rule['log_dir'],fnamebase+'.jobs')
    payload = {'rule':rule,'shell_list':shell_list}

    with open(jobfile,'w') as f:
        json.dump(payload,f,sort_keys=False,indent=2)
        f.write('\n')

    return fnamebase

def execute_command(config,job_numb,cmd,env):
    'execute command locally'
    fname = f'{config["log_prefix"]}{timestamp_now()}.{config["name"]}'
    foutname = os.path.join(config['log_dir'],fname+'.out')
    ferrname = os.path.join(config['log_dir'],fname+'.err')
    fout = open(foutname,'w')
    ferr = open(ferrname,'w')

    failed = False
    print(f'executing job number {job_numb+1} locally:',end='')
    sys.stdout.flush()
    try:
        subprocess.run(cmd,env=env,shell=True,check=True,stdout=fout,stderr=ferr)
    except subprocess.CalledProcessError:
        failed = True

    fout.close()
    ferr.close()

    if failed: print(' failed')
    else:      print(' okay')

    return failed

# def write_qsub_file(rule,qsub_script,jobname,njobs,jobfile):
#     with open(qsub_script,'w') as f:
#         f.write("#!/bin/bash -l\n")
#         f.write(f"#$ -N {jobname}\n")
#         f.write(f"#$ -cwd\n")
#         f.write(f"#$ -V\n")
#         f.write(f"#$ -t 1-{njobs}\n")
#         f.write(f"#$ -o {rule['log_dir']}/{jobname}.$TASK_ID.out\n")
#         f.write(f"#$ -e {rule['log_dir']}/{jobname}.$TASK_ID.err\n")
#         f.write(f"#$ -l h_rt={rule['time']}\n")
#         f.write(f"#$ -l mem={rule['mem']}\n")
#         f.write(f"#$ -l tmpfs={rule['tmpfs']}\n")
#         f.write(f"#$ -pe {rule['pe']} {rule['cores']}\n")

#         #run payload through fakemake
#         f.write(f"{sys.executable} {sys.argv[0]} --qsub {jobfile}\n")

def write_qsub_file(rule,qsub_script,jobname,njobs,jobfile):
    'fill out the qsub job script template and write to file ready to pass to qsub'

    if rule['qsub_template'] == 'default':
        f_in = open(os.path.join(os.path.dirname(__file__),"qsub_template.sh"))
    else:
        f_in = open(rule['qsub_template'])

    f_out = open(qsub_script,'w')

    env = copy.deepcopy(rule)
    env["njobs"] = njobs
    env["jobname"] = jobname
    env["jobfile"] = jobfile
    env["python"] = sys.executable
    env["fakemake"] = sys.argv[0]

    for line in f_in: f_out.write(line.format(**env))

    f_out.close()
    f_in.close()

def submit_job_qsub(rule,shell_list,job_list):
    '''
    issue the qsub command to spawn an array of jobs
    wait for completion before checking the status of each job
    '''
    #rule['python_executable'] = sys.executable
    #rule['python_path'] = sys.path
    jobname = write_jobfile(rule,shell_list)
    qsub_script = os.path.join(rule['log_dir'],jobname+'.qsub')
    jobfile = os.path.join(rule['log_dir'],jobname+'.jobs')
    njobs = len(shell_list)

    write_qsub_file(rule,qsub_script,jobname,njobs,jobfile)

    cmd = f"qsub -sync y {qsub_script}"

    env = copy.deepcopy(os.environ)

    foutname = os.path.join(rule['log_dir'],jobname+'.out')
    ferrname = os.path.join(rule['log_dir'],jobname+'.err')
    fout = open(foutname,'w')
    ferr = open(ferrname,'w')

    failed = False
    print(f'executing qsub job array with {njobs} jobs:',end='')
    sys.stdout.flush()
    try:
        subprocess.run(cmd,env=env,shell=True,check=True,stdout=fout,stderr=ferr)
    except subprocess.CalledProcessError:
        failed = True

    fout.close()
    ferr.close()

    #delay to allow for shared filesystem latency on status files
    time.sleep(int(rule['remote_delay_secs']))

    if failed: print(' failed')
    else:      print(' okay')

    #check each individual job's status file
    for job_numb,item in enumerate(shell_list):
        status_file = f'{rule["log_dir"]}/{jobname}.{job_numb+1}.status'

        job_failed = False
        if not os.path.exists(status_file):
            print(f"job {job_numb+1} probably failed to start")
            job_failed = True
        else:
            with open(status_file) as f: status = f.read().strip()
            print(f"job {job_numb+1} status: {status}")
            if status != "okay": job_failed= True

        if not job_failed: continue

        #deal with output of failed job
        handle_failed_outputs(rule,job_list[job_numb]["output"])

def qsub_execute_job(jobfile):
    '''
    runs on a compute node under grid engine control
    execute one job of a job array spawned by qsub
    invoked automatically by fakemake using the --qsub option
    '''
    rule,shell_list = read_jobfile(jobfile)

    assert jobfile.endswith('.jobs')

    status_file = f'{jobfile[:-5]}.{os.environ["SGE_TASK_ID"]}.status'
    job_numb = int(os.environ["SGE_TASK_ID"]) - 1 #convert from 1 to 0 based
    item = shell_list[job_numb]

    cmd = generate_full_command(rule,item)
    env = generate_job_environment(rule,job_numb,len(shell_list))

    #created/truncate to record that the job is starting
    f = open(status_file,'w')
    f.close()

    failed = False
    try:
        subprocess.run(cmd,env=env,shell=True,check=True)
    except subprocess.CalledProcessError:
        failed = True

    #record outcome of the job
    f = open(status_file,'w')
    if failed:
        f.write('failed\n')
        f.close()
        sys.exit(1)
    else:
        f.write('okay\n')
        f.close()
        sys.exit(0)

def execute_jobs(rule,shell_list,job_list):
    if rule['exec'] == 'local':
        #local serial execution
        for job_numb,item in enumerate(shell_list):
            cmd = generate_full_command(rule,item)
            env = generate_job_environment(rule,job_numb,len(shell_list))
            failed = execute_command(rule,job_numb,cmd,env)

            if failed: handle_failed_outputs(rule,job_list[job_numb]["output"])

    elif rule['exec'] == 'qsub':
        #qsub execution using an array job
        submit_job_qsub(rule,shell_list,job_list)

    else:
        raise Exception(f"unsupported execution method {rule['exec']}")

def timestamp_now():
    return datetime.datetime.fromtimestamp(time.time()).strftime("%Y%m%d.%H%M%S.%f")

def read_jobfile(fname):
    'load the rule config and shell_list from json'

    if not fname.endswith('.jobs'):
        fname += '.jobs'

    with open(fname) as f:
         payload = json.load(f)

    return payload['rule'], payload['shell_list']

def process_rule(config,rule):
    #validate rule and substitute placeholders
    input,output,shell = setup_rule(config,rule)

    print(f'\nrule: {rule["name"]}')

    #fill in any globbing placeholders in inputs and outputs
    #return list of values, one per potential job
    job_list = generate_job_list(rule,input,output)

    #determine if all inputs are present
    #and if outputs need to be regenerated
    job_list,shell_list = generate_shell_commands(rule,job_list,shell)

    print(f'{len(shell_list)} jobs generated')

    if len(shell_list) == 0: return

    #check current working directory agrees with configured value
    check_cwd(rule)

    #execute jobs locally or remotely from the jobfile
    execute_jobs(rule,shell_list,job_list)
