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
import itertools

global_state = {}

#regex patterns to match non nested {%placeholders} and {$environment variables}
environ_regx  = r'\{\$.+?\}'       #{$name}
scalar_regx   = r'\{%.+?\}'        #{%scalar} or {%list[]}
glob2job_regx  = r'\{\*.+?\}'      #{*glob} ==> split into separate jobs
glob2list_regx  = r'\{\+.+?\}'     #{+glob} ==> become in-job list
glob2both_regx  = r'\{[\+\*].+?\}' #{+/*glob} ==> either type of glob placeholder
list2job_regx  = r'\{=.+?\}'       #{=list} ==> split into separate jobs
list2list_regx  = r'\{-.+?\}'      #{-list} ==> become in-job list

warning_prefix = 'WARNING: '

default_global_config =\
{
    'fm':{
        'prefix':'fm-',                  #log file/job name  prefix: qsub doesn't like numerical names
        'log_dir':'fakemake_logs',       #name of subfolder for logging
        'remote_delay_secs':'10',        #wait this long after remote jobs incase of latency
        'stale_output_file':'ignore',    #ignore,delete,recycle (also applies to symlinks)
        'stale_output_dir':'ignore',     #ignore,delete,recycle
        'failed_output_file':'stale',    #delete,recycle,stale,ignore (also applies to symlinks)
        'failed_output_dir':'stale',     #delete,recycle,stale,ignore    
        'missing_parent_dir':'create',   #ignore,create
        'recycle_bin':'recycle_bin',     #name of recycle bin folder
        'job_count':'FM_NJOBS',          #env variable: how many jobs spawned by current action
        'job_number':'FM_JOB_NUMBER',    #env variable: 1 based job numbering within the current action
        'bash_setup':'source ~/.bashrc\nset -euo pipefail\nset +o history',
        'conda_setup':'',
    },
    'qsub':{
        'template':'default',        #template job script: "default" or path to your own
        'time':'02:00:00',           #$ -l h_rt={time}
        'mem':'4G',                  #$ -l mem={mem}
        'tmpfs':'10G',               #$ -l tmpfs={tmpfs}
        'pe':'smp',                  #$ -pe {pe} {cores}
        'cores':'1',                 #$ -pe {pe} {cores}
        'log_dir':'{%fm/log_dir}'    #log dir for qsub stdout and stderr files
                                     #$ -o/e {log_dir}/{jobname}.$TASK_ID.out/err
    },
    'exec':'local',                  #default execution environment
}

class Conf:
    def __init__(self,src=None):
        self.d = None

        if src == None:
            #initialise to empty
            self.d = {}

        elif src == "defaults":
            #initialise to defaults
            self.d = copy.deepcopy(default_global_config)

        elif type(src) == dict:
            #copy another dict
            self.d = copy.deepcopy(src)
            self.stringify()

        elif type(src) == Conf:
            #copy another Conf
            self.d = copy.deepcopy(src.getdict())

        else:
            raise Exception(f"unsupported source class {type(src)}")

    def __contains__(self,key):
        try:
            self[key]
            return True
        except:
            return False

    def stringify(self):
        #make sure any numerical types are converted to string
        for key,value in self.items():
            if type(value) != str:
                self[key] = str(value)

    def includes_and_loads(self,parent_file):
        'action any includes or file loads'

        parent_path = os.path.dirname(parent_file)

        while True:
            #find next remaining include, process it, then start again
            #to avoid changing dict while iterating it
            load_type = None

            for subkeys,value in self.subkey_items():
                if 'includes' in subkeys:
                    load_type = 'include'

                    if not subkeys[-2] == 'includes':
                        raise Exception(f'misplaced "includes" directive at {"/".join(subkeys)}')

                    if not value.endswith('.yml'):
                        raise Exception(f'all "includes" must be .yml files for {"/".join(subkeys)}')

                    filename = value

                elif value.startswith('{>'):
                    load_type = 'file'
                    assert value.endswith('}'),f'misformed file load placeholder {value} at {"/".join(subkeys)}'

                    filename = value[2:-1]

                #exit loop as soon as we find an include or file load
                if load_type != None: break

            #got to the end without finding anything therefore we are done
            if load_type == None: return

            #action the include or load
            if load_type == 'include':
                self.do_include(subkeys,filename,parent_path)
            else:
                self.do_load(subkeys,filename,parent_path)
    
    def do_load(self,subkeys,filename,parent_path):
        'load a list of values from the file'
        key = '/'.join(subkeys)

        #extract any parameters
        sep = None
        row = None
        col = None

        fname = filename        
        if '[' in filename:
            assert filename[-1] == ']',f'misformed load parameters in {filename} at {key}'
            params = filename.split('[')[-1][:-1]
            fname = filename.split('[')[0]

            if 'R' in params:
                tok = params.split('R')
                assert len(tok) == 2,f'misformed load parameters in {filename} at {key}'
                sep = tok[0]
                row = int(tok[1])

            elif 'C' in params:
                tok = params.split('C')
                assert len(tok) == 2,f'misformed load parameters in {filename} at {key}'
                sep = tok[0]
                col = int(tok[1])

            else:
                raise Exception(f'misformed load parameters in {filename} at {key}')

        #get full path to the included file
        full_path = make_fullpath(fname,parent_path)

        if not os.path.exists(full_path):
            raise Exception(f"missing file {full_path} at {key}")

        data = []

        counter = -1
        with open(full_path) as f:
            if sep == None:
                data = f.read().strip()
            else:
                for line in f:
                    counter += 1

                    if row != None and row == counter:
                        data = line.strip().split(sep)
                        break

                    if col != None:
                        values = line.strip().split(sep)
                        data.append(values[col])

        assert data != [],f'failed to extract any data from {filename}'

        self[subkeys] = data

    def do_include(self,subkeys,filename,parent_path):
        #where to put the loaded variables
        base_keys = subkeys[:-2]

        #delete the include item so it won't be actioned a second time
        self.delete(subkeys)

        #delete the includes list if now empty
        if len(self.getitem(subkeys[:-1])) == 0:
            self.delete(subkeys[:-1])

        self.do_include_inner(filename,parent_path,base_keys)

    def do_include_inner(self,filename,parent_path,base_keys):
        #get full path to the included file
        full_path = make_fullpath(filename,parent_path)

        #load the YAML
        with open(full_path) as f:
            yml = yaml.safe_load(f)
            assert type(yml) == dict, f'include {full_path} must be a dict at top level'
            self.update(yml,base=base_keys)

    def update(self,src,base=None):
        '''
        merge in all items from another Conf object or dict
        values of src take priority (overwrite) any duplicate keys in self
        in lists shared between self and src, the list in self gets emptied
        before the items from src are added
        place the new items at the specified base key
        default is the root key
        '''

        if type(src) == dict: src = Conf(src)

        if type(src) != Conf:
            raise Exception(f"unsupported object type {type(src)}")

        if base:
            base_key,base_subkeys = self.get_key_and_subkeys(base)

        for src_key,src_value in src.containers():
            if base:
                self_key = '/'.join((base_key,src_key))
            else:
                self_key = src_key

            try:
                self_value = self[self_key]
            except:
                self_value = None

            if self_value == None:
                if type(src_value) == list:
                    #pad out to required new length
                    self[self_key] = [None] * len(src_value)
                elif type(src_value) == dict:
                    self[self_key] = {}
                else:
                    raise Exception(f"unsupported container type {type(src_value)} at {src_key}")

            elif type(src_value) == list:
                #wipe any existing list values, pad out to required new length
                self[self_key] = [None] * len(src_value)
            elif type(src_value) == dict and type(self_value) != dict:
                #only create empty dict if not already a dict
                #therefore preserver any existing keys
                self[self_key] = {}

        for src_key,src_value in src.items():
            if base:
                self_key = '/'.join((base_key,src_key))
            else:
                self_key = src_key

            self[self_key] = src_value

    def override(self,src):
        '''
        merge in all items from another Conf object or dict
        values of self take priority for any duplicate keys
        '''

        tmp = Conf(src)
        tmp.update(self)
        self.update(tmp)

    def getdict(self):
        return self.d

    def show(self,label=None,indent=2):
        if label: print(label)
        print(json.dumps(self.d,sort_keys=False,indent=indent))
        print()

    def keys(self,item=None):
        'generator function to crawl through all leaf keys'

        for key,value in self.items(): yield key

    def items(self,item=None):
        '''
        generator function to crawl through all leaf key:value pairs
        return key as '/'.join(subkey_list)
        rather than a list of subkeys
        '''

        for subkeys,value in self.subkey_items(item):
            yield '/'.join(subkeys),value

    def containers(self,item=None):
        '''
        generator function to crawl through all key:value pairs
        returning only containers not leaves
        return key as '/'.join(subkey_list)
        rather than a list of subkeys
        '''

        for subkeys,value in self.subkey_containers(item):
            yield '/'.join(subkeys),value

    def subkey_containers(self,item=None):
        '''
        generator function to crawl through all key:value pairs
        returning only containers not leaves
        return key as a list of subkeys
        rather than '/'.join(subkey_list)
        '''

        if item == None: item = self.d

        if type(item) == dict:
            for key in item:
                if type(item[key]) not in (list,dict): continue

                yield [key],item[key]

                for subkey,value in self.subkey_containers(item[key]):
                    yield [key]+subkey,value

        elif type(item) == list:
            for i in range(len(item)):
                if type(item[i]) not in (list,dict): continue

                yield [str(i)],item[i]

                for subkey,value in self.subkey_containers(item[i]):
                    yield [str(i)]+subkey,value

        else:
            raise Exception(f"cannot generate items from type {type(item)}")

    def subkey_items(self,item=None):
        '''
        generator function to crawl through all leaf key:value pairs
        return key as explicit list of subkeys
        '''

        if item == None: item = self.d

        if type(item) == dict:
            for key in item:
                if type(item[key]) not in (list,dict):
                    yield [key],item[key]
                else:
                    for subkey,value in self.subkey_items(item[key]):
                        yield [key]+subkey,value
        elif type(item) == list:
            for i in range(len(item)):
                if type(item[i]) not in (list,dict):
                    yield [str(i)],item[i]
                else:
                    for subkey,value in self.subkey_items(item[i]):
                        yield [str(i)]+subkey,value
        else:
            raise Exception(f"cannot generate items from type {type(item)}")

    def get_key_and_subkeys(self,key):
        '''
        given either type of key
        return the '/'.join and list forms of the key
        '''

        if type(key) == list:
            for x in key: assert type(x) == str
            key = '/'.join(key)
            subkeys = key.split('/')
        elif type(key) == str:
            subkeys = key.split('/')
        else:
            raise Exception(f'invalid key type {type(key)}')

        assert len(subkeys) > 0, f"invalid empty key {key}"

        return key,subkeys

    def __getitem__(self,key):
        '''
        lookup an item by hierachical key
        error if key does not exist
        '''

        return self.getitem(key,delete=False)

    def delete(self,key):
        'delete the item pointed to by the key'
        self.getitem(key,delete=True)

    def getitem(self,key,delete=False):
        '''
        lookup an item by hierachical key
        key can be as 'subkey1/subkey2...'
        or ['subkey1','subkey2'...]
        for lists the "subkey" is the normal python index
        or N for list length
        or any single non digit for joined version of whole list
        for dicts the "subkey" is a normal key string
        ??or '' for list of key strings
        error if key does not exist
        delete=True means delete instead of returning value
        '''

        key,subkeys = self.get_key_and_subkeys(key)

        item = self.d
        final = False
        for i,sk in enumerate(subkeys):
            if i == len(subkeys)-1: final = True

            if type(item) == dict:
                if sk == '':
                    #special meaning: get list of keys
                    assert not delete, "cannot delete keys using the empty key"
                    if final:
                        return list(item.keys())
                    else:
                        item = list(item.keys())
                        continue

                assert sk in item, f"key error for {sk} within {key}"

                if not final:
                    item = item[sk]
                    continue

                if delete:
                    del item[sk]
                    return None

                return item[sk]

            elif type(item) == list:
                if final and not is_valid_int(sk):
                    #special subscript meaning
                    assert not delete,f'invalid use of special subscript {sk} at {key}'
                    return self.special_subscript(sk,item)

                try:
                    item[int(sk)]
                except:
                    raise Exception(f"key error for {sk} within {key}")

                sk = int(sk)

                if not final:
                    item = item[sk]
                    continue

                if delete:
                    del item[sk]
                    return None

                return item[sk]
            else:
                raise Exception(f"unsupported data type in {key}")

    def special_subscript(self,sk,item):
        if sk == 'N':
            return str(len(item))
        elif sk == '\\t':
            return '\t'.join(item)
        elif sk == '\\n':
            return '\n'.join(item)
        else:
            return sk.join(item)

    def __setitem__(self,key,value):
        '''
        set item by hierachical key
        creating any parent nodes required
        '''

        key,subkeys = self.get_key_and_subkeys(key)

        item = self.d
        assign = False
        for i,sk in enumerate(subkeys):
            if i == len(subkeys) - 1: #leaf reached, time to assign the new value
                assign = True

            if type(item) == dict:
                assert sk != '', "cannot assign to dict key list"

                if assign:
                    item[sk] = value
                    return

                assert sk in item, "this implementation doesn't support implicit parent node creation"

                #trying to avoid using the below as it prevents ints being
                #used as dict keys - since the below assumes they are list indexes
                # if not sk in item:
                #     #need to create an empty container
                #     try:
                #         int(subkeys[i+1]) #see if next subkey implies a list or a dict
                #         item[sk] = []
                #     except:
                #         item[sk] = {}

                item = item[sk]
            elif type(item) == list:
                try:
                    #validate the list index
                    sk = int(sk)
                except:
                    raise Exception(f"key error for {sk} within {key}")

                #assigning single new item in existing list
                assert sk < len(item) and sk >= -len(item)

                if assign:
                    item[sk] = value
                    return

                item = item[sk]

            else:
                raise Exception(f"unsupported data type in {key}")

    def sub_vars(self,src=None):
        '''
        substitute all {$...} and {%...} placeholders
        or throw exception if 10 iterations is not enough
        '''

        #where to get the replacement values from
        if src == None: src = self

        counter = 10
        while True:
            changed = False
            counter -= 1
            if self.sub_values(os.environ,environ_regx): changed = True
            if self.sub_values(src,scalar_regx):         changed = True
            if not changed: break
            assert counter > 0, 'unable to resolve all placeholders in 10 iterations'

    def sub_values(self,src,regx):
        '''
        substitute all placeholders
        eg {$environ}, {%scalars}
        '''

        changed = False

        for key,value in self.items():
            new_val,flag = self.do_subs(src,regx,value)
            if flag:
                self[key] = new_val
                changed = True

        return changed

    def do_subs(self,src,regx,value):
        '''
        sub all simple placeholders in value from src
        src can be a dict (os.environ) or a Conf 
        '''

        changed = False

        while True:
            m = re.search(regx,value)
            if m == None: break #no more matches

            key = m.group(0)[2:-1]   #{+key}, {$key}, {%hierachical/key}, {%list/0}
            sub = src[key]

            if type(sub) != str:
                raise Exception(f'key {key} resolves to a {type(sub)} not a string')

            value = value[:m.start(0)] + sub + value[m.end(0):]

            changed = True

        return value,changed

    def make_log_dir(self):
        if not os.path.exists(self['fm/log_dir']):
            message(f'creating missing log_dir {self["fm/log_dir"]}')

            #note: still creating missing log_dir in dryrun mode
            #so we can record messages and create qsub scripts
            #for inspection by the user
            #therefore not using the makedirs wrapper function
            os.makedirs(self['fm/log_dir'])

def makedirs(path):
    if is_active(): os.makedirs(path)

def make_fullpath(filename,parent_path):
    'return full path to the include file path'

    if filename.startswith('./'):
        #wrt current working directory
        full_path = filename
    else:
        #wrt parent_path (probably the main yaml file)
        full_path = os.path.join(parent_path,filename)

    return full_path

def is_valid_int(subscript):
    try:
        int(subscript)
        return True
    except:
        return False

def valid_subscript(subscript,list_var):
    if not is_valid_int(subscript):
        return False

    if int(subscript) < -len(list_var) or int(subscript) >= len(list_var):
        return False

    return True

def show(item,label=None,indent=2):
    if label: print(label)
    print(json.dumps(item,sort_keys=False,indent=indent))
    print()

def check_cwd(config):
    if 'working_dir' in config:
        if os.path.realpath(os.getcwd()) != config['working_dir']:
            print("warning: current path is not the expected working directory")
            print(f"expecting: {config['working_dir']}")
            print(f"but found {os.path.realpath(os.getcwd())}")

def toplevel_include(config,src,path):
    '''
    include into toplevel of config from a single path string
    or list of path strings
    '''

    if type(src) == str:
        config.do_include_inner(src,path,base_keys=[])

    elif type(src) == list:
        for src_path in src:
            assert type(src_path) == str, "invalid include list, list items must be strings"
            config.do_include_inner(src_path,path,base_keys=[])

    else:
        raise Exception('invalid include item, must be string or list of strings')

    config.update(item[item_type])
    #config.show("loaded config")
    config.includes_and_loads(path)
    #config.show("actioned includes and loads")

def is_active():
    '''
    return True if activity_state() is "active"
    return False if activity_state() is "dryrun" or "inactive"
    '''

    if activity_state() == 'active': return True
    return False

def activity_state():
    '''
    return the current activity state: active, dryrun or inactive
    active = run everything as normal (action is active, dryrun is off)
    dryrun = pretend to run (action is active, but dryrun is on)
    inactive = run nothing (action is not active)
    '''

    if global_state['is_active'] == False:
        #do not run the current action at all
        return 'inactive'

    if global_state['args'].dry_run == True:
        #run the current action in dryrin mode only
        return 'dryrun'

    #run as normal
    return 'active'

def update_activity_state(action):
    '''
    update the global_state['is_active'] state based on the action name
    this implements the run-only,run-from and run-until command line options
    '''

    #incase the config has changed since the last action
    if 'fm/log_dir' in action: global_state['log_dir'] = action['fm/log_dir']

    #see if we've reached the run-from rule
    if global_state['args'].run_from and action['name'] == global_state['args'].run_from:
        global_state['reached_run_from'] = True

    #see if we've reached the run-until rule
    if global_state['args'].run_until and action['name'] == global_state['args'].run_until:
        global_state['reached_run_until'] = True

    #if run-only is active then the action must be in the approved list
    if global_state['args'].run_only and not action['name'] in global_state['args'].run_only:
        global_state['is_active'] = False

    #see if we're within the run-from to run-until interval
    elif global_state['reached_run_from'] and not global_state['reached_run_until']:
        global_state['is_active'] = True
    else:
        global_state['is_active'] = False

def init_global_state(args,config):
    global global_state

    #only set global state once
    if global_state != {}: return

    assert args != None

    #store dry_run, run_only, run_from, run_until settings
    global_state['args'] = args

    #store changeable states
    #run-from implies we start off not running any action yet
    if args.run_from:
        global_state['reached_run_from'] = False
    else:
        global_state['reached_run_from'] = True

    #becomes true when we encounter the named run-until action, if any
    global_state['reached_run_until'] = False

    #used to generate the log file names
    global_state['start_time'] = timestamp_now()

    global_state['is_active'] = None
    global_state['log_dir'] = config['fm/log_dir']

    update_activity_state({'name':None})


def load_pipeline(path,parent_file):
    'path is relative to the parent path'

    parent_path = os.path.dirname(parent_file)
    full_path = make_fullpath(path,parent_path)

    assert full_path != parent_file
    return parse_yaml(full_path),full_path

def find_duplicates(conf_list):
    all_keys = set()

    for config in conf_list:
        for key in config.keys():
            assert not key in all_keys, f"duplicate key {key}"
            all_keys.add(key)

def fill_out_all(config,action,path):
    'filling in missing placeholders etc'
    action.override(config)
    action.includes_and_loads(path)
    action.sub_vars()
    action.make_log_dir()

def setup_action(config,action,path):
    #separate out and canonicalise input, output and shell
    action,input,output,shell = split_action(action)

    #merge config into action but action has priority
    fill_out_all(config,action,path)

    #try to ensure placeholder subtitution is not ambiguous
    #find_duplicates([action,input,output])

    #substitute all placeholders except for shell
    #which likely contains input/output variables yet to be determined
    input.sub_vars(src=action)
    output.sub_vars(src=action)

    return action,input,output,shell

def split_action(action):
    '''
    separate out input, output and shell from action
    convert all to Conf
    note:action begins as a simple dict
    '''

    input = action['input']
    output = action['output']
    shell = action['shell']

    del action['input']
    del action['output']
    del action['shell']

    #convert simple form input/output into dictionary form
    if type(input) == str:  input = {'input':input}
    if type(output) == str: output = {'output':output}
    shell = { 'shell':shell }

    action = Conf(action)
    input = Conf(input)
    output = Conf(output)
    shell = Conf(shell)

    return action,input,output,shell

def parse_yaml(fname):
    'parse yaml into list of items'
    with open(fname) as f:
        result = yaml.safe_load(f)

    assert type(result) == list
    
    return result

def generate_job_list(action,input,output):
    '''
    expand list and glob placeholders in input patterns to final file names
    this also expands the job list to one job per matching combination
    fill out output patterns with the matched names
    '''

    #seed job "list" containing just the unexpanded input(s) and output(s)
    #which may have placeholders in need of expansion
    job_list = [{"input":input,"output":output}]

    #expand in-job list patterns
    for job in job_list: generate_list2list(action,job)

    # print("after injob list")
    # for job in job_list: show_job(job)

    #expand separate-jobs list patterns
    new_list = []
    for job in job_list: new_list += generate_list2jobs(action,job)
    job_list = new_list

    if len(job_list) == 0: return []

    # print("after between list")
    # for job in job_list: show_job(job)

    #expanded input pattern globs (both the injob and separate jobs kind)
    new_list = []
    for job in job_list: new_list += generate_glob_jobs(action,job)
    job_list = new_list

    if len(job_list) == 0: return []

    # print("after globs")
    # for job in job_list: show_job(job)

    return job_list

def expand_lists(action,item,key,value):
    #get list of injob list expansion placeholders
    ph_names = {}
    for m in re.finditer(list2list_regx,value):
        name = m.group(0)[2:-1]   #"{-name}" ==> "name"
        ph_names[name] = True

    #if no injob list placeholders present then this key does
    #not need to be expanded into a new sublist
    if(len(ph_names)) == 0: return

    #expand patterns into lists of patterns
    #filled out with the values from the list(s)
    #where multiple lists are present expand to all combinations of the lists
    meta_list = [ action[key] for key in ph_names ]

    new_list = []

    #"product" implements nested for loops over all the lists
    for value_list in itertools.product(*meta_list):
        src = { name:value_list[i] for i,name in enumerate(ph_names.keys()) }

        dst = Conf({"value":value})
        #substitute in the correct values from the list variables
        dst.sub_values(src,list2list_regx)
        new_list.append(dst["value"])

    item[key] = new_list

def generate_list2list(action,job):
    '''
    expand in-job list placeholders for each job inplace
    no new jobs are created, only lists within jobs are expanded
    if they contain {-name} type placeholder(s)
    '''

    #identify all in-job type placeholders {-name} in input/output patterns

    for key,value in job["input"].items():
        expand_lists(action,job["input"],key,value)

    for key,value in job["output"].items():
        expand_lists(action,job["output"],key,value)

def generate_list2jobs(action,job):
    '''
    expand between-job list placeholders to one independent job
    per placeholder value for single lists
    or value combination where multiple placeholders are in the same pattern
    '''

    #identify all listjob type placeholders {=name} in input/output patterns
    ph_names = {}

    for key,value in job["input"].items():
        for m in re.finditer(list2job_regx,value):
            name = m.group(0)[2:-1]   #"{=name}" ==> "name"
            if name not in ph_names: ph_names[name] = True

    for key,value in job["output"].items():
        for m in re.finditer(list2job_regx,value):
            name = m.group(0)[2:-1]   #"{=name}" ==> "name"
            if name not in ph_names: ph_names[name] = True

    #nothing to expand if no list patterns present in input patterns
    if len(ph_names) == 0:
        job["job_vars"] = Conf()
        return [ job ]

    #expand input and output patterns into complete paths
    #(excepting any remaining glob placeholders)
    #filled out with the values from the list(s)
    #where multiple lists are present expand to all combinations of the lists
    new_list = []
    meta_list = [action[key] for key in ph_names]

    #do nested for loops over all the lists
    for value_list in itertools.product(*meta_list):
        src = { name:value_list[i] for i,name in enumerate(ph_names.keys()) }

        new_input = Conf(job["input"])
        new_output = Conf(job["output"])

        #substitute in the correct values from the list variables
        new_input.sub_values(src,list2job_regx)
        new_output.sub_values(src,list2job_regx)

        #add new matches to job
        new_lists = Conf(src)
        new_list.append({"input":new_input,"output":new_output,"job_vars":new_lists})

    return new_list

def build_glob_and_regx(pattern):
    input_glob = ''
    input_regx = '^'
    ph_dict = {}
    prev_end = 0
    names1 = {} # in-job names
    names2 = {} # between-job names

    #convert glob-type placeholders into proper glob and regex query formats
    for m in re.finditer(glob2both_regx,pattern):
        name = m.group(0)[2:-1]   #{+name} or {*name}
        start = m.start(0)

        input_glob += pattern[prev_end:start] + '*'
        input_regx += re.escape(pattern[prev_end:start])

        #{*name} defines a set of separate jobs
        if name not in ph_dict:
            ph_dict[name] = True
            input_regx += '(?P<' + name + '>[^/]+)'
        else:
            #back reference to previous job or list placeholder
            input_regx += '(?P=' + name + ')'

        if m.group(0)[1] == '+':
            names1[name] = True
        else:
            names2[name] = True

        prev_end = m.end(0)

    input_glob += pattern[prev_end:]
    input_regx += re.escape(pattern[prev_end:]) + '$'
    
    return input_glob,input_regx,names1,names2

def expand_into_lists(item):
    '''
    find all keys which contain {+name} type placeholders
    expand these into lists
    unless they are already members of a list
    return a dict with all the {+name} keys in
    set to True if it was expanded
    False otherwise
    '''

    key_dict = {}
    for subkeys,pattern in item.subkey_items():
        #print('/'.join(subkeys))

        if re.search(glob2list_regx,pattern):
            key = '/'.join(subkeys)

            if len(subkeys) == 1 or type(item[subkeys[:-1]]) != list:
                #expand into list
                key_dict[key] = True
            else:
                #already a list, do not expand further
                key_dict[key] = False

    for key in key_dict:
        if key_dict[key] == True:
            item[key] = [ item[key] ]

    return key_dict

def make_alt_key(all_between):
    #store the new job keyed by the list of "between" values
    #ie *not* the normal input file key
    #as we may need to append new injob list items to already created jobs
    #and need to look up by the list of between values
    for x in all_between.values(): assert not '/' in x
    sorted_keys = list(all_between.keys())
    sorted_keys.sort()
    job_key = '/'.join([all_between[key] for key in sorted_keys])
    #note: job_key will be '' if no between keys present
    #but it's still a usable dict key

    return job_key

def generate_glob_jobs(action,job):
    '''
    glob filename to match the two glob type placeholders
    {+name} are globs that generate in-job lists, eg:
    info_file: "samples/{+name}.txt"
    becomes a list within a job like
    info_file:
      - "samples/frog.txt"
      - "samples/toad.txt" etc
    {*name} are globs that generate separate jobs ("varies between jobs"), eg:
    info_file: "samples/{*name}.txt"
    becomes separate jobs where:
    info_file: "samples/frog.txt"
    and then
    info_file: "samples/toad.txt" in the next job etc
    '''

    # print("before expand into lists")
    # show_job(job)
    #exit()

    #expand into lists all input and output keys which contain injob list placeholders
    expanding_input_keys = expand_into_lists(job["input"])
    expanding_output_keys = expand_into_lists(job["output"])

    # print("after expand into lists")
    # show_job(job)
    #exit()

    #glob matches keyed by input pattern name
    glob_matches = {}

    #glob each input file pattern separately
    for key,pattern in job["input"].items():
        #convert the pattern-with-placeholders into a glob and regex string
        #names1 contains all the injob glob keys {+name}
        #names2 contains all the between job glob keys {*name}
        input_glob,input_regx,names1,names2 = build_glob_and_regx(pattern)

        # print(input_glob,input_regx)
        # show(names1)
        # show(names2)
        #exit()

        glob_matches[key] = []
        
        #find paths using iglob, match to placeholders using regex
        for path in glob.iglob(input_glob):
            m = re.fullmatch(input_regx,path)
            assert m is not None,"regex cannot extract placeholders from the globbed path!"

            #split the re matches into the two types of placeholder
            injob = {}
            between = {}
            for name in m.groupdict():
                if name in names1:
                    injob[name] = m.groupdict()[name] #{+name}
                elif name in names2:
                    between[name] = m.groupdict()[name] #{*name}
                else:
                    raise Exception(f"missing name {name}")

            #store the match as the completed path plus the placeholder values
            glob_matches[key].append([key,path,injob,between])
        
        #no jobs are generated if any input pattern has zero matches
        if len(glob_matches[key]) == 0: return []

        #sort into alphabetical order by path
        glob_matches[key].sort(key=lambda item: item[1])

    # print("glob_matches")
    # for key in glob_matches:
    #     show(glob_matches[key])
    #exit()

    #find glob_matches where the placeholders
    #have matching values across all input patterns (including injob and between)
    #reject all others
    new_job_dict = {}
    meta_list = [ glob_matches[key] for key in glob_matches ]

    # show(meta_list)
    # exit()

    #do nested for loops over all the lists in meta_list
    #to generate all-vs-all combinations
    for value_list in itertools.product(*meta_list):
        #build the potential new jobs
        all_injob = {}  #accumulate all injob variables across all input patterns
        all_between = {}#accumulate all between job variables across all input patterns
        #inputs = {}
        new_input = Conf(job["input"])

        conflict = False
        for key,path,injob,between in value_list:
            # print(key)
            # print(path)
            # print("injob")
            # show(injob)
            # print("between")
            # show(between)
            # print()

            #find conflicts for any variables shared between paths
            if find_conflicts(all_injob,injob) or find_conflicts(all_between,between):
                #print("conflict")
                conflict = True
                break

            # store completed path under the input file key
            # note key is in '/'.join(subkeys) format
            new_input[key] = path 

        #mismatching combo of input path variables therefore no job generated
        if conflict:
            #print('conflict')
            continue

        #all_injob and all_between now contain the final placeholder values
        #with conflicts removed

        # print(job_key)
        # exit()

        #input placeholders are already filled out by the glob
        ##new_input = Conf()
        ##new_input.update(inputs)

        #fill in placeholders in output patterns here
        new_output = Conf(job['output'])
        for key,pattern in new_output.items():
            new_output.sub_values(all_between,glob2job_regx)
            new_output.sub_values(all_injob,glob2list_regx)

        #make alternative job key based on "between"
        job_key = make_alt_key(all_between)
        appending = False
        if not job_key in new_job_dict:
            #create new job upon first encounter with this particular combo of "between" values
            new_job = \
            {
                "input":new_input,
                "output":new_output,
                "job_vars":Conf(job["job_vars"]),
                "glob_vars":Conf(all_between),
            }
            new_job["glob_lists"] = Conf({k:[v] for k,v in injob.items()})
            new_job_dict[job_key] = new_job

        else:
            #we're adding new injob list members to an existing job
            appending = True
            new_job = new_job_dict[job_key]
            for exp_key in expanding_input_keys:
                # print(exp_key)
                # new_job['input'].show()
                if type(new_job['input'][exp_key]) == list:
                    # print(exp_key)
                    new_job['input'][exp_key].append(new_input[exp_key][0])
                else:
                    parent_key = '/'.join(exp_key.split('/')[:-1])
                    # print('>',parent_key,'<')
                    new_job['input'][parent_key].append(new_input[exp_key])
                # new_job['input'].show()

            for exp_key in expanding_output_keys:
                # print(exp_key)
                # new_job['output'].show()

                if type(new_job['output'][exp_key]) == list:
                    # print(exp_key)
                    new_job['output'][exp_key].append(new_output[exp_key][0])
                else:
                    parent_key = '/'.join(exp_key.split('/')[:-1])
                    # print('>',parent_key,'<')
                    new_job['output'][parent_key].append(new_output[exp_key])
                # new_job['output'].show()

            for exp_key in injob:
                new_job['glob_lists'][exp_key].append(injob[exp_key])

            # exit()

        # show_job(new_job_dict[job_key])
        # exit()

    #convert new_jobs into list
    job_list = [new_job_dict[k] for k in new_job_dict]

    #for job in job_list: show_job(job)

    return job_list

def show_job(job):
    print('=====job=====')
    for key in job:
        print(key)
        job[key].show()
    print()

def find_conflicts(all_vals,values):
    'return True upon first conflict found'
    for name in values:
        if name not in all_vals:
            all_vals[name] = values[name]
        elif all_vals[name] != values[name]:
            return True #conflict found
    return False #no conflicts found

def check_input_mtimes(input):
    'return newest mtime if all inputs present otherwise return None'

    #verify all input paths present
    missing = False
    newest_mtime = 0
    for item,path in input.items():
        if not os.path.exists(path):
            message(f'missing input {path}')
            missing = True
            continue

        if is_stale(path):
            message(f'stale input {path}')
            missing = True
            continue

        mtime = os.path.getmtime(path)
        if mtime > newest_mtime: newest_mtime = mtime

    if missing: return None

    return newest_mtime

def check_output_mtimes(output):
    'return oldest mtime if all outputs present otherwise return None'

    #dummy future time guaranteed not to be older than a real file
    oldest_mtime = time.time() + 31e6 

    for item,path in output.items():
        if not os.path.exists(path):
            return None

        mtime = os.path.getmtime(path)
        if mtime < oldest_mtime: oldest_mtime = mtime

    
    return oldest_mtime

def generate_shell_commands(action,job_list,shell):

    for job_numb,job in enumerate(job_list):
        #check all inputs present, return newest mtime
        newest_input = check_input_mtimes(job['input'])

        #one or more inputs missing, job not runnable
        if newest_input == None:
            #flag job for removal from the list
            job_list[job_numb] = None
            warning(f'one or more inputs missing or stale, job not runable')
            continue

        #all inputs present
        #check if any outputs missing, older than newest input
        #more marked as stale
        oldest_output = check_output_mtimes(job['output'])

        #all outputs present and newer than newest input: no need to run
        #unless run: always is set for this action 
        if oldest_output != None and oldest_output > newest_input:
            if 'run' in action and action['run'] == 'always':
                message('"run always" option forcing outputs to be treated as stale')
            else:
                #flag job for removal from the list
                job_list[job_numb] = None
                continue

        #remove stale output files/symlinks/directories before action is run
        handle_stale_outputs(action,job['output'])

        #create any missing output *parent* directories
        if action['fm/missing_parent_dir'] == 'create':
            create_output_dirs(action,job['output'])

    #remove non-runable jobs
    job_list = [job for job in job_list if job is not None]
    shell_list = []

    for job_numb,job in enumerate(job_list):
        #merge input and output filenames and placeholders into action variables
        config = Conf(action)
        config.update(job['input'])
        config.update(job['output'])

        #config.show()

        #substitute remaining placeholders in shell command
        shell_final = Conf(shell)
        shell_final.sub_values(os.environ,environ_regx)
        shell_final.sub_values(config,scalar_regx)
        shell_final.sub_values(job['job_vars'],list2job_regx)
        shell_final.sub_values(job['glob_vars'],glob2job_regx)
        shell_final.sub_values(job['glob_lists'],glob2list_regx)
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
    
    recycle_bin = config['fm/recycle_bin']

    assert path != recycle_bin

    new_name = os.path.join(recycle_bin,os.path.basename(path))

    #make sure not to overwrite existing recycle bin file
    if os.path.exists(new_name):
        new_name += ".%d"%random.randrange(1000000)
        if os.path.exists(new_name):
            raise Exception(f"File {new_name} already in recycle bin")

    #move to recycle bin
    move_item(path,new_name)

def move_item(path,new_name):
    if is_active():
        shutil.move(path,new_name)

def create_output_dirs(config,outputs):
    '''
    if the parent folder of any output is missing create it
    '''

    prev_parent = None

    for item in outputs.keys():
        path = outputs[item]
        parent = os.path.dirname(path)
        if os.path.exists(parent): continue

        #avoid message spamming about missing directories in dry-run mode
        if parent != prev_parent:
            message(f'creating missing output directory {parent}')

        makedirs(parent)
        prev_parent = parent

def remove_item(path):
    if is_active(): os.remove(path)

def remove_tree(path):
    if is_active(): shutil.rmtree(path)

def warning(item,end='\n'):
    message(warning_prefix+item,end=end)

def message(item,end='\n'):
    item = timestamp_now_nice() + ' ' + str(item)

    if global_state['args'].quiet != True:
        print(item,end=end)
        sys.stdout.flush()

    if global_state['args'].nologs != True:
        path = os.path.join(global_state['log_dir'],global_state['start_time']+'.messages')

        try:
            with open(path,'a') as f:
                f.write(item+end)
                f.close()
        except:
            print(timestamp_now_nice() + ' ' + warning_prefix+f'log file {path} not writeable')

def handle_stale_outputs(config,outputs):
    '''
    the job is going to be rerun therefore any outputs already present
    are treated as stale and recycled or deleted
    '''

    for item in outputs.keys():
        path = outputs[item]
        if not os.path.exists(path): continue

        if os.path.islink(path) or os.path.isfile(path):
            if config['fm/stale_output_file'] == 'delete':
                message(f'deleting stale output file {path}')
                remove_item(path)
            elif config['fm/stale_output_file'] == 'recycle':
                message(f'recycling stale output file {path}')
                recycle_item(config,path)
            elif config['fm/stale_output_file'] == 'ignore':
                message(f'ignoring stale output file {path}')
                continue
            else:
                raise Exception(f'unknown option for stale_output_file: {config["fm/stale_output_file"]}')

        elif os.path.isdir(path):
            if config['fm/stale_output_dir'] == 'delete':
                remove_tree(path)
            elif config['fm/stale_output_dir'] == 'recycle':
                recycle_item(config,path)
            elif config['fm/stale_output_dir'] == 'ignore':
                continue
            else:
                raise Exception(f'unknown option for stale_output_dir: {config["fm/stale_output_dir"]}')
                
        else:
            raise Exception(f'unsupported output type {path}')

def handle_failed_outputs(config,outputs):
    '''
    the job has failed therefore any outputs present are untrustworthy
    and should be removed, recycled, flagged as old etc
    '''

    for item in outputs.keys():
        path = outputs[item]
        if not os.path.exists(path): continue

        if os.path.islink(path) or os.path.isfile(path):
            if config['fm/failed_output_file'] == 'delete':
                remove_item(path)
            elif config['fm/failed_output_file'] == 'recycle':
                recycle_item(config,path)
            elif config['fm/failed_output_file'] == 'stale':
                make_stale(path)
            elif config['fm/failed_output_file'] == 'ignore':
                continue
            else:
                raise Exception(f'unknown option for failed_output_file: {config["fm/failed_output_file"]}')

        elif os.path.isdir(path):
            if config['fm/failed_output_dir'] == 'delete':
                remove_tree(path)
            elif config['fm/failed_output_dir'] == 'recycle':
                recycle_item(config,path)
            elif config['fm/failed_output_dir'] == 'stale':
                make_stale(path)
            elif config['fm/failed_output_dir'] == 'ignore':
                continue
            else:
                raise Exception(f'unknown option for failed_output_dir: {config["fm/failed_output_dir"]}')
                
        else:
            raise Exception(f'unsupported output type {path}')

def generate_full_command(config,shell):
    'generate the bash commands to run before the job commands'

    cmd_list = []

    if 'fm/bash_prefix' in config and config['fm/bash_prefix'] != '':
        cmd_list.append(config['fm/bash_prefix'])

    if 'conda' in config and config['conda'] != '':
        if 'fm/conda_setup_command' in config and config['fm/conda_setup_command'] != '':
            cmd_list.append(config['fm/conda_setup_command'])
        cmd_list.append(f'conda activate {config["conda"]}')

    cmd_list.append(shell)

    return '\n'.join(cmd_list)

def generate_job_environment(config,job_numb,njobs):
    'copy local environment with a few adjustments'

    env = copy.deepcopy(os.environ)
    env[ config['fm/job_count'] ] = f'{njobs}'
    env[ config['fm/job_number'] ] = f'{job_numb+1}' # convert from 0 to 1 based to match SGE_TASK_ID

    return env

def write_jobfile(action,shell_list):
    'save action config and shell_list as json'

    fnamebase = f'{action["fm/prefix"]}{timestamp_now()}.{action["name"]}'
    jobfile = os.path.join(action['fm/log_dir'],fnamebase+'.jobs')
    payload = {'action':action.getdict(),'shell_list':shell_list}

    with open(jobfile,'w') as f:
        json.dump(payload,f,sort_keys=False,indent=2)
        f.write('\n')

    return fnamebase

def execute_command(config,job_numb,cmd,env):
    'execute command locally'
    fname = f'{config["fm/prefix"]}{timestamp_now()}.{config["name"]}'
    foutname = os.path.join(config['fm/log_dir'],fname+'.out')
    ferrname = os.path.join(config['fm/log_dir'],fname+'.err')

    message(f'job {job_numb+1} executing locally...')

    if cmd.endswith('\n'): message(f'{cmd}',end='')
    else:                  message(f'{cmd}')

    if not is_active():
        #dry-run: signal job completed ok without running it
        message(f'skipping execution')
        return False

    failed = False
    fout = open(foutname,'w')
    ferr = open(ferrname,'w')

    try:
        subprocess.run(cmd,env=env,shell=True,check=True,stdout=fout,stderr=ferr)
    except subprocess.CalledProcessError:
        failed = True

    fout.close()
    ferr.close()

    if failed: warning('command failed')
    else:      message('command completed normally')

    return failed

def write_qsub_file(action,qsub_script,jobname,njobs,jobfile):
    'fill out the qsub job script template and write to file ready to pass to qsub'

    if action['qsub/template'] == 'default':
        f_in = open(os.path.join(os.path.dirname(__file__),"qsub_template.sh"))
    else:
        f_in = open(action['qsub/template'])

    if is_active():
        message(f'creating qsub jobscript {qsub_script}')
    else:
        message(f'creating qsub jobscript {qsub_script} (even in dry-run mode)')

    f_out = open(qsub_script,'w')

    env = copy.deepcopy(action['qsub'])
    env["njobs"] = njobs
    env["jobname"] = jobname
    env["jobfile"] = jobfile
    env["python"] = sys.executable
    env["fakemake"] = sys.argv[0]

    for line in f_in: f_out.write(line.format(**env))

    f_out.close()
    f_in.close()

def submit_job_qsub(action,shell_list,job_list):
    '''
    issue the qsub command to spawn an array of jobs
    wait for completion before checking the status of each job
    '''

    jobname = write_jobfile(action,shell_list)
    qsub_script = os.path.join(action['fm/log_dir'],jobname+'.qsub')
    jobfile = os.path.join(action['fm/log_dir'],jobname+'.jobs')
    njobs = len(shell_list)

    write_qsub_file(action,qsub_script,jobname,njobs,jobfile)

    cmd = f"qsub -sync y {qsub_script}"

    env = copy.deepcopy(os.environ)

    message(f'executing {cmd} with {njobs}...')

    something_failed = False

    if is_active():
        foutname = os.path.join(action['fm/log_dir'],jobname+'.out')
        ferrname = os.path.join(action['fm/log_dir'],jobname+'.err')
        fout = open(foutname,'w')
        ferr = open(ferrname,'w')
        #sys.stdout.flush()

        try:
            subprocess.run(cmd,env=env,shell=True,check=True,stdout=fout,stderr=ferr)
        except subprocess.CalledProcessError:
            something_failed = True

        fout.close()
        ferr.close()

        if something_failed: warning(f'qsub call failed')
        else:                message(f'qsub call succeeded')

    #delay to allow for shared filesystem latency on status files
    message(f'sleeping for {action["fm/remote_delay_secs"]} seconds to allow for filesystem latency...')
    if activity_state() == 'active': time.sleep(int(action['fm/remote_delay_secs']))

    #check each individual job's status file and expected outputs
    #some jobs may not have failed even if qsub returned an error code?
    for job_numb,item in enumerate(shell_list):
        status_file = f'{action["fm/log_dir"]}/{jobname}.{job_numb+1}.status'

        failed = False
        if not os.path.exists(status_file):
            if activity_state() == 'dryrun':
                warning(f"job {job_numb+1} produced no status file {status_file} due to dryrun mode")
            else:
                warning(f"job {job_numb+1} produced no status file {status_file} and probably failed to start")
            failed = True
        else:
            with open(status_file) as f: status = f.read().strip()

            if status == "okay":
                message(f"job {job_numb+1} status: {status}")
            else:
                warning(f"job {job_numb+1} status: {status}")
                failed = True

        if failed == False:
            #job has also failed if not all required outputs were created
            failed = verify_expected_outputs(action,job_list[job_numb]["output"],job_list[job_numb]["input"])

        if failed == True:
            #deal with output of failed job
            handle_failed_outputs(action,job_list[job_numb]["output"])
            something_failed = True

    return something_failed

def qsub_execute_job(jobfile):
    '''
    runs on a compute node under grid engine control
    execute one job of a job array spawned by qsub
    invoked automatically by fakemake using the --qsub option
    '''
    action,shell_list = read_jobfile(jobfile)

    action = Conf(action)

    #action.show()

    assert jobfile.endswith('.jobs')

    spawn_job(jobfile,action,shell_list)

def spawn_job(jobfile,action,shell_list):
    'spawn user job'

    status_file = f'{jobfile[:-5]}.{os.environ["SGE_TASK_ID"]}.status'
    job_numb = int(os.environ["SGE_TASK_ID"]) - 1 #convert from 1 to 0 based
    item = shell_list[job_numb]

    cmd = generate_full_command(action,item)
    env = generate_job_environment(action,job_numb,len(shell_list))

    #create/truncate to record that the job is starting
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

def verify_expected_outputs(action,outputs,inputs):
    '''
    verify that all specified output files exist
    and are newer than the job start time
    if any are missing or stale flag job as failed

    implies jobs must use touch to update any output file that does not need altering
    or else exclude it from the list of required outputs
    '''

    #find newest input file
    newest_mtime = check_input_mtimes(inputs)

    #if any input(s) missing cannot verify that outputs are not stale
    #implies jobs cannot delete their inputs
    #ie this would need to be handled by a separate job
    if newest_mtime == None:
        warning('unable to verify all outputs are fresh due to missing input(s)')
        return True

    failed = False
    for item in outputs.keys():
        path = outputs[item]
        if not os.path.exists(path):
            warning(f'missing output {path}')
            failed = True 
        elif is_stale(path) or os.path.getmtime(path) < newest_mtime:
            message(f'stale output {path}')
            failed = True

    if failed:
        #job has failed
        return True

    #signal job seems to have completed okay
    message('all outputs look good')
    return False

def execute_jobs(action,shell_list,job_list):
    something_failed = False

    if action['exec'] == 'local':
        #local serial execution
        for job_numb,item in enumerate(shell_list):
            cmd = generate_full_command(action,item)
            env = generate_job_environment(action,job_numb,len(shell_list))

            #reports only explicit job failure from exit code
            failed = execute_command(action,job_numb,cmd,env)

            if failed == False:
                #job has also failed if not all required outputs were created
                failed = verify_expected_outputs(action,job_list[job_numb]["output"],job_list[job_numb]["input"])

            if failed == True:
                #carry out config controlled action on outputs of failed job
                #ie delete, recycle, mark as stale or ignore
                handle_failed_outputs(action,job_list[job_numb]["output"])
                something_failed = True

    elif action['exec'] == 'qsub':
        #qsub execution using an array job
        something_failed = submit_job_qsub(action,shell_list,job_list)

    else:
        raise Exception(f"unsupported execution method {action['exec']}")

    return something_failed

def timestamp_now():
    return datetime.datetime.fromtimestamp(time.time()).strftime("%Y%m%d.%H%M%S.%f")

def timestamp_now_nice():
    return datetime.datetime.fromtimestamp(time.time()).strftime("%Y-%m-%d %H:%M:%S")

def read_jobfile(fname):
    'load the action config and shell_list from json'

    if not fname.endswith('.jobs'):
        fname += '.jobs'

    with open(fname) as f:
         payload = json.load(f)

    return payload['action'], payload['shell_list']

def process_action(config,action,path):
    #validate action and substitute placeholders
    action,input,output,shell = setup_action(config,action,path)
    #action.show("setup action")

    #see if the rule name triggers a change
    update_activity_state(action)

    if activity_state() == 'inactive':
        message(f'skipping action {action["name"]}')
        return False

    if 'run' in action and action['run'] == 'never':
        message(f'skipping action {action["name"]}')
        return False

    if activity_state() == 'dryrun':
        message(f'dry-running action {action["name"]}')
    else:
        message(f'running action: {action["name"]}')

    #generate list of potential jobs by filling out glob and list placeholders
    #this step does not pay any attention to file time stamps
    job_list = generate_job_list(action,input,output)
    message(f'placeholder expansion found {len(job_list)} potential jobs')
    if len(job_list) == 0: return False

    #determine if all inputs are present and non-stale
    #determine if any outputs need (re)generating
    job_list,shell_list = generate_shell_commands(action,job_list,shell)
    message(f'input/output file checking found {len(shell_list)} runable jobs')
    if len(shell_list) == 0: return False

    #check current working directory agrees with configured value
    check_cwd(action)

    #execute jobs locally or remotely from the jobfile
    if execute_jobs(action,shell_list,job_list):
        warning(f'action {action["name"]} had failed job(s)')
        #signal action failed
        return True

    #signal action completed ok
    return False

def process(pipeline,path,config=None,args=None):
    '''
    process a YAML pipeline from the file args.yaml
    '''

    #initially set config to default values if none provided
    if config == None:
        #set conf to defaults
        config = Conf(src="defaults")

    #on first call store args and init stateful variables that
    #implement run-only, run-from and run-until options
    init_global_state(args,config)

    counter = 0

    message(f"starting pipeline {path}")

    while counter < len(pipeline):
        item = pipeline[counter]
        counter += 1
        assert type(item) == dict
        assert len(item) == 1
        item_type = list(item.keys())[0]

        if item_type == 'action':
            #show(item[item_type],item_type)
            if process_action(config,item[item_type],path):
                if activity_state() == 'dryrun':
                    warning('action failed due to dryrun mode, continuing pipeline anyway')
                else:
                    warning('aborting pipeline due to failed action')
                    return 2

        elif item_type == 'config':
            #add new config to the existing one, overriding any shared keys
            #config.show("orig config")
            message(f'updating config')
            config.update(item[item_type])
            #config.show("loaded config")
            config.includes_and_loads(path)
            #config.show("actioned includes and loads")

        elif item_type == 'include':
            #toplevel include: load and insert the yaml items in place of the include item
            message(f'including {item[item_type]}')
            new_pipeline,new_path = load_pipeline(item[item_type],path) #new_path ignored
            counter -= 1
            del pipeline[counter]
            pipeline = pipeline[:counter] + new_pipeline + pipeline[counter:]

        elif item_type == 'module':
            #process a nested pipeline without affecting the config of any
            #following items
            new_pipeline,new_path = load_pipeline(item[item_type],path)
            sub_config = Conf(config)
            process(new_pipeline,new_path,config=sub_config)

        else:
            raise Exception(f"unsupported item type {item_type}")

    message(f"reached end of pipeline {path}")

    return 0
