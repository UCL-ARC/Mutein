import yaml
import re
import copy
import json
import os
import glob
import datetime
import time
import shutil

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
    'temp_job_dir':'fakemake_jobs',
    'exec':'local',
    'conda':'fm_default_env',
    'stale_output_file':'ignore',   #ignore,delete,recycle also applies to symlinks
    'stale_output_dir':'ignore',    #ignore,delete,recycle
    'missing_parent_dir':'create',     #ignore,create
    'recycle_bin':'recycle_bin',
    'array_size_variable':'FM_ARRAY_SIZE', #how many jobs in current array
    #'normalize_paths':'true',               #true,false
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

def process(config,pipeline):
    check_cwd(config)

    #pipeline yaml must be a list of items all called "rule"
    for i,item in enumerate(pipeline):
        assert type(item) == dict
        assert len(item) == 1
        item_type = list(item.keys())[0]
        assert item_type == 'rule'
        #show(item[item_type],item_type)
        process_rule(config,item[item_type])

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

def setup_config(pipeline):
    config = copy.deepcopy(default_global_config)

    #include the config from the pipeline yaml file if present
    item = pipeline[0]
    assert type(item) == dict and len(item) == 1
    item_type = list(item.keys())[0]
    if item_type == 'config': config.update(item[item_type])

    #check all keys and values are simple strings
    #check for forbidden keys
    for key,value in config.items():
        assert type(key) == str
        assert type(value) == str
        assert key not in nonconfig_keys

    #substitute any environment variables
    sub_environ(config)

    #substitute any fakemake variables
    sub_vars(config)

    show(config,"config")

    del pipeline[0]

    if not os.path.exists(config['temp_job_dir']):
        os.makedirs(config['temp_job_dir'])

    return config

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

    #print(primary)
    #print(input_glob)
    #print(input_regx)

    job_list = []

    #find paths using iglob, match to placeholders using regex
    for path in glob.iglob(input_glob):
        #print(path)
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
        if not os.path.exists(path):
            print(f"missing input {path}")
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
            print(f"missing output {path}")
            all_outputs_present = False
            continue

        mtime = os.path.getmtime(path)
        if mtime < oldest_mtime: oldest_mtime = mtime

    #return oldest mtime if all outputs present otherwise None
    if all_outputs_present == False: return None
    return oldest_mtime

def generate_final_shell_commands(rule,job_list,shell):
    name = rule['name']
    for job_numb,job in enumerate(job_list):
        #check all inputs present, returned newest mtime
        newest_input = check_input_mtimes(job['input'])

        #one or more inputs missing, job not runnable
        if newest_input == None:
            #flag job for removal from the list
            job_list[job_numb] = None
            continue

        #all inputs present
        print(f"rule:{name} job:{job_numb} all inputs present")
        newest_datetime = datetime.datetime.fromtimestamp(newest_input)
        print(f"newest input {newest_datetime}")

        #check if any outputs missing or older than newest input
        oldest_output = check_output_mtimes(job['output'])

        #all outputs present and newer than newest input: no need to run
        if oldest_output != None and oldest_output > newest_input:
            print(f"rule:{name} job:{job_numb} all outputs present")
            oldest_datetime = datetime.datetime.fromtimestamp(oldest_output)
            print(f"oldest output {oldest_datetime}")

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
        show(shell_final,"shell")
        shell_list.append(shell_final["shell"])

    return shell_list

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
        print(path)
        parent = os.path.dirname(path)
        if os.path.exists(parent): continue
        print(f"creating {parent}")
        os.makedirs(parent)

def handle_stale_outputs(config,outputs):
    '''
    the job is going to be rerun therefore any output already present
    is treated as stale and recycled or deleted
    '''

    print("handle stale outputs")

    for item in outputs:
        path = outputs[item]
        print(path)
        if not os.path.exists(path): continue

        if os.path.islink(path) or os.path.isfile(path):
            if config['stale_output_file'] == 'delete':
                print(f"remove file {path}")
                os.remove(path)
            elif config['stale_output_file'] == 'recycle':
                print(f"recycle file {path}")
                recycle_item(config,path)
            elif config['stale_output_file'] == 'ignore':
                print(f"ignore file {path}")
                continue
            else:
                raise Exception(f'unknown option for stale_output_file: {config["stale_output_file"]}')

        elif os.path.isdir(path):
            if config['stale_output_dir'] == 'delete':
                print(f"remove dir {path}")
                shutil.rmtree(path)
            elif config['stale_output_dir'] == 'recycle':
                print(f"recycle dir {path}")
                recycle_item(config,path)
            elif config['stale_output_dir'] == 'ignore':
                print(f"ignore dir {path}")
                continue
            else:
                raise Exception(f'unknown option for stale_output_dir: {config["stale_output_dir"]}')
                
        else:
            raise Exception(f'unsupported output type {path}')

def execute_jobs(rule,jobfile):
    with open(jobfile) as f:
        shell_list = json.load(f)

    if rule['exec'] == 'local':
        for item in shell_list:
            print("local exec:\n",item)
        pass
        #local execution
    elif rule['exec'] == 'qsub':
        pass
        #qsub execution
    else:
        raise Exception(f"unsupported execution method {rule['exec']}")

def write_joblist(rule,shell_list):
    timestamp = datetime.datetime.fromtimestamp(time.time()).strftime(".%Y%m%d-%H%M%S.%f")
    jobfile = os.path.join(rule['temp_job_dir'],rule['name']+timestamp)

    with open(jobfile,'w') as f:
        json.dump(shell_list,f,sort_keys=False,indent=2)
        f.write('\n')

    return jobfile

def process_rule(config,rule):
    #validate rule and substitute placeholders
    input,output,shell = setup_rule(config,rule)

    print('------------------------')
    show(rule,"rule")
    show(input,"input")
    show(output,"output")
    show(shell,"shell")

    #fill in any globbing placeholders in inputs and outputs
    #return list of values, one per job
    job_list = generate_job_list(rule,input,output)
    show(job_list,"jobs")

    if len(job_list) == 0:
        print(f"rule {rule['name']} no jobs generated")
        return

    #determine if all inputs are present
    #and if outputs need to be regenerated
    shell_list = generate_final_shell_commands(rule,job_list,shell)

    #write list of shell commands to jobfile
    jobfile = write_joblist(rule,shell_list)

    #execute jobs locally or remotely
    execute_jobs(rule,jobfile)
