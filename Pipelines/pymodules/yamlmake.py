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

default_config_str =\
'''
    ym:
        prefix:             'ym-'           #log file/job name  prefix: qsub doesn't like numerical names
        remote_delay_secs:  '10'            #wait this long after remote jobs incase of latency
        stale_output_file:  'ignore'        #ignore,delete,recycle (also applies to symlinks)
        stale_output_dir:   'ignore'        #ignore,delete,recycle
        failed_output_file: 'stale'         #delete,recycle,stale,ignore (also applies to symlinks)
        failed_output_dir:  'stale'         #delete,recycle,stale,ignore    
        missing_parent_dir: 'create'        #ignore,create
        recycle_bin:        'recycle_bin'   #name of recycle bin folder
        job_count:          'YM_NJOBS'      #env variable: how many jobs spawned by current action
        job_number:         'YM_JOB_NUMBER' #env variable: 1 based job numbering within the current action
        conda_setup:        ''              #run just before trying to activate the conda env the command requested
        conda_prefix:       ''              #a prefix to apply to the name of every conda environment
        #run before every shell action
        bash_setup: |
          source ~/.bashrc
          set -euo pipefail
          set +o history 
    qsub:
        template:           'default'       #template job script: "default" or path to your own
        time:               '02:00:00'      #$ -l h_rt={time}
        mem:                '4G'            #$ -l mem={mem}
        tmpfs:              '10G'           #$ -l tmpfs={tmpfs}
        pe:                 'smp'           #$ -pe {pe} {cores}
        cores:              '1'             #$ -pe {pe} {cores}
'''

default_global_config = yaml.safe_load(default_config_str)

#global meta configuration
meta = {}

#regex patterns to match non nested {%placeholders} and {$environment variables}
environ_regx  = r'\{\$.+?\}'       #{$name}
scalar_regx   = r'\{%.+?\}'        #{%scalar} or {%list[]}
glob2job_regx  = r'\{\*.+?\}'      #{*glob} ==> split into separate jobs
glob2list_regx  = r'\{\+.+?\}'     #{+glob} ==> become in-job list
glob2both_regx  = r'\{[\+\*].+?\}' #{+/*glob} ==> either type of glob placeholder
list2job_regx  = r'\{=.+?\}'       #{=list} ==> split into separate jobs
list2list_regx  = r'\{-.+?\}'      #{-list} ==> become in-job list

warning_prefix = 'WARNING: '
error_prefix = 'ERROR: '
illegal_chrs = '\\\n\r\t\'" *?:;,/#%&{}<>+`|=$!@'

default_log_dir = 'yamlmake_logs'

col=\
{
    'none':'\033[0m',
    'black':'\033[1;30m',
    'red':'\033[1;31m',
    'green':'\033[1;32m',
    'yellow':'\033[1;33m',
    'blue':'\033[1;34m',
    'purple':'\033[1;35m',
    'cyan':'\033[1;36m',
    'white':'\033[1;37m',
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

    def is_empty(self):
        'return True if contains no keys'

        if self.d == {}: return True
        return False

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
        #optional extra info can follow the filename in square brackets
        #if not square brackets entire file is loaded including multiple lines
        if '[' in filename:
            #filename[<options>]
            #options are: <sep>R<numb> or <sep>C<numb>
            #sep is the field separator
            #R = extract a row
            #C = extract a column
            #numb = which row or column to extract
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

    def no_placeholders(self):
        'verify no placeholders remain'

        reg_list = [
            environ_regx,
            scalar_regx,
            glob2job_regx,
            glob2list_regx,
            list2job_regx,
            list2list_regx,
        ]

        for key,value in self.items():
            for regx in reg_list:
                if re.search(regx,value):
                    raise Exception(f'unmatched placeholder {value} found in {key}')

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

    # log_dir now either default or set on command line
    # otherwise cannot be changed
    # def make_log_dir(self):
    #     if not os.path.exists(self['ym/log_dir']):
    #         message(f'creating missing log_dir {self["ym/log_dir"]}')

    #         #note: still creating missing log_dir in dryrun mode
    #         #so we can record messages and create qsub scripts
    #         #for inspection by the user
    #         #therefore not using the makedirs wrapper function
    #         os.makedirs(self['ym/log_dir'])

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
        if os.path.realpath(os.getcwd()) != os.path.realpath(config['working_dir']):
            print("warning: current path is not the expected working directory")
            print(f"expecting: {config['working_dir']}")
            print(f"but found {os.path.realpath(os.getcwd())}")

# def toplevel_include(config,src,path):
#     '''
#     include into toplevel of config from a single path string
#     or list of path strings
#     '''

#     if type(src) == str:
#         config.do_include_inner(src,path,base_keys=[])

#     elif type(src) == list:
#         for src_path in src:
#             assert type(src_path) == str, "invalid include list, list items must be strings"
#             config.do_include_inner(src_path,path,base_keys=[])

#     else:
#         raise Exception('invalid include item, must be string or list of strings')

#     config.update(item[item_type])
#     #config.show("loaded config")
#     config.includes_and_loads(path)
#     #config.show("actioned includes and loads")

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

    if meta['is_active'] == False:
        #do not run the current action at all
        return 'inactive'

    if meta['args'].dry_run == True:
        #run the current action in dryrin mode only
        return 'dryrun'

    #run as normal
    return 'active'

def update_activity_state(action):
    '''
    update the meta['is_active'] state based on the action name
    this implements the run-only,run-from and run-until command line options
    '''

    #incase the config has changed since the last action
    if 'ym/prefix' in action:
        meta['prefix'] = action['ym/prefix']

    #see if we've reached the run-from rule
    if meta['args'].run_from and action['name'] == meta['args'].run_from:
        meta['reached_run_from'] = True

    #see if we've reached the run-until rule
    if meta['args'].run_until and action['name'] == meta['args'].run_until:
        meta['reached_run_until'] = True

    #if run-only is active then the action must be in the approved list
    if meta['args'].run_only and not action['name'] in meta['args'].run_only:
        meta['is_active'] = False

    #see if we're within the run-from to run-until interval
    elif meta['reached_run_from'] and not meta['reached_run_until']:
        meta['is_active'] = True
    else:
        meta['is_active'] = False

def init_meta(args,config):
    global meta

    #only set global state once
    if meta != {}: return

    assert args != None

    #store dry_run, run_only, run_from, run_until settings
    meta['args'] = args

    #store changeable states
    #run-from implies we start off not running any action yet
    if args.run_from:
        meta['reached_run_from'] = False
    else:
        meta['reached_run_from'] = True

    #becomes true when we encounter the named run-until action, if any
    meta['reached_run_until'] = False

    #used to generate the log file names
    meta['start_time'] = timestamp_now()

    meta['is_active'] = None
    meta['prefix'] = config['ym/prefix']

    if args.log_dir != None:
        meta['log_dir'] = args.log_dir
    else:
        meta['log_dir'] = default_log_dir

    #note: still creating missing log_dir in dryrun mode
    #so we can record messages and create qsub scripts
    #for inspection by the user
    #therefore not using the makedirs wrapper function
    if not os.path.exists(meta['log_dir']):
        os.makedirs(meta['log_dir'])
        message(f"created missing log_dir {meta['log_dir']}")

    #update with a fake action name that will not match any rule
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
    '''
    merge action with global config
    do any includes and loads
    filling in missing placeholders
    create log folder if missing
    '''

    action.override(config)
    action.includes_and_loads(path)
    action.sub_vars()
    #action.make_log_dir()

def validate_action(action):
    'check for required fields'

    assert 'name' in action, 'all actions must have a name field'
    assert 'shell' in action, 'all actions must have a shell field'

    #assert 'output' in action, 'all actions must have an output field'
    #assert 'input' in action, 'all actions must have an input field'
    
    if not 'input' in action: action['input'] = {}   #allow empty input field
    if not 'output' in action: action['output'] = {} #allow empty output field

    #name is used to general file and job names
    #therefore prevent any characters that cause problems for filenames
    for ch in illegal_chrs:
        assert not ch in action['name'], f'cannot use "{ch}" in action name, consider a description field for longer text'

    assert action['name'][0] != '.', f'action name {action["name"]} cannot start with a dot'

def setup_action(config,action,path):
    #validate the action
    validate_action(action)

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

    #decided to remove this little convenience as it makes
    #explaining input/output more complex
    # #convert simple form input/output into dictionary form
    # if type(input) == str:  input = {'input':input}
    # if type(output) == str: output = {'output':output}

    shell = { 'shell':shell }
    action = Conf(action)
    input = Conf(input)
    output = Conf(output)
    shell = Conf(shell)

    return action,input,output,shell

def parse_yaml(fname):
    'parse yaml into list of items'

    #make sure all shell: fields are literal block scalars
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("shell:"): continue
            if not line == "shell: |":
                raise Exception('all shell fields must use the literal block scalar style "shell: |" to preserve newlines')

    with open(fname) as f:
        result = yaml.safe_load(f)

    #allow empty files, such as when all items are commented out
    if result == None: result = []

    assert type(result) == list, f"YAML file {fname} does not contain a list of YAMLmake items"
    
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

    #expand separate-jobs list patterns
    new_list = []
    for job in job_list: new_list += generate_list2jobs(action,job)
    job_list = new_list

    if len(job_list) == 0: return []

    #expanded input pattern globs (both the injob and separate jobs kind)
    new_list = []
    for job in job_list: new_list += generate_glob_jobs(action,job)
    job_list = new_list

    if len(job_list) == 0: return []

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
    expand job spawning list placeholders to one independent job
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
    expand these into lists unless they are already members of a list
    return a dict with all the {+name} keys in
    set to True if it was expanded
    False otherwise
    '''

    #find what needs expanding
    key_dict = {}
    for subkeys,pattern in item.subkey_items():
        if re.search(glob2list_regx,pattern):
            key = '/'.join(subkeys)

            if len(subkeys) == 1 or type(item[subkeys[:-1]]) != list:
                #expand into list
                key_dict[key] = True
            else:
                #already a list, do not expand further
                key_dict[key] = False

    #perform the expansion
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

def generate_glob_matches(job_input):
    #signal that no input patterns were specified
    if job_input.is_empty(): return "no-input"

    glob_matches = {}

    #glob each input field file pattern separately
    for key,pattern in job_input.items():
        #convert the pattern-with-placeholders into a glob and regex string
        #names1 contains all the injob glob keys {+name}
        #names2 contains all the between job glob keys {*name}
        input_glob,input_regx,names1,names2 = build_glob_and_regx(pattern)

        glob_matches[key] = []
        
        #find paths using iglob, match to placeholders using regex
        for path in glob.iglob(input_glob):
            m = re.fullmatch(input_regx,path)
            
            if m == None:
                #regex cannot extract placeholders
                #usually indicates two *'s in the glob matches different strings
                #where as in the regex those two are required to be the same
                #eg {*accession}/{*accession}
                # ==> glob is */*
                # ==> one path generated is A123/A123.txt
                # which doesn't match the regex as "A123" != "A123.txt"
                continue

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
        #this includes input patterns without placeholders where the
        #hardcoded path/file is missing
        if len(glob_matches[key]) == 0:
            message(f'no match for input {key}: {pattern}')
            return "no-match"

        #sort into alphabetical order by path
        glob_matches[key].sort(key=lambda item: item[1])

    return glob_matches

def generate_glob_jobs(action,job):
    '''
    glob filenames to match the two glob type placeholders
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

    #expand into lists all input and output keys which contain injob list placeholders
    #ready to receive the expanded items later
    expanding_input_keys = expand_into_lists(job["input"])
    expanding_output_keys = expand_into_lists(job["output"])

    #find glob matches keyed by input pattern name
    #returns "no-match" if no matches found to at least one input pattern
    #return "empty" if no inputs patterns specified
    glob_matches = generate_glob_matches(job["input"])

    #no jobs generated from this input pattern
    if glob_matches == "no-match": return []

    if glob_matches == "no-input":
        #no glob expansion to perform, return job as is with empty glob variable dicts
        new_job = \
        {
            "input":Conf(job['input']),
            "output":Conf(job['output']),
            "job_vars":Conf(job['job_vars']),
            "glob_vars":Conf(),
            "glob_lists":Conf(),
        }
        new_job['input'].no_placeholders()
        new_job['output'].no_placeholders()
        return [new_job]

    #find glob_matches where the placeholders
    #have matching values across all input patterns (including injob and between)
    #reject all others
    new_job_dict = {}
    meta_list = [ glob_matches[key] for key in glob_matches ]

    #do nested for loops over all the lists in meta_list
    #to generate all-vs-all combinations
    for value_list in itertools.product(*meta_list):
        #build the potential new jobs
        all_injob = {}  #accumulate all injob variables across all input patterns
        all_between = {}#accumulate all between job variables across all input patterns
        new_input = Conf(job["input"])

        conflict = False
        for key,path,injob,between in value_list:

            #find conflicts for any variables shared between paths
            if find_conflicts(all_injob,injob) or find_conflicts(all_between,between):
                conflict = True
                break

            # store completed path under the input file key
            # note key is in '/'.join(subkeys) format
            new_input[key] = path 

        #mismatching combo of input path variables therefore no job generated
        if conflict:
            continue

        #all_injob and all_between now contain the final placeholder values
        #with conflicts removed

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
                if type(new_job['input'][exp_key]) == list:
                    new_job['input'][exp_key].append(new_input[exp_key][0])
                else:
                    parent_key = '/'.join(exp_key.split('/')[:-1])
                    new_job['input'][parent_key].append(new_input[exp_key])

            for exp_key in expanding_output_keys:
                if type(new_job['output'][exp_key]) == list:
                    new_job['output'][exp_key].append(new_output[exp_key][0])
                else:
                    parent_key = '/'.join(exp_key.split('/')[:-1])
                    new_job['output'][parent_key].append(new_output[exp_key])

            for exp_key in injob:
                new_job['glob_lists'][exp_key].append(injob[exp_key])

    #convert new_jobs into list
    job_list = [new_job_dict[k] for k in new_job_dict]

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
    '''
    return newest mtime if all inputs present
    otherwise return "missing" if any missing
    or "empty" if there are no inputs listed
    '''

    #no inputs to check
    if input.is_empty(): return "empty"

    #verify all input paths present
    missing = False
    newest_mtime = None

    for item,path in input.items():
        if not os.path.exists(path):
            #should get filtered out by generate_glob_jobs already
            message(f'missing input {path}')
            missing = True
            continue

        if is_stale(path):
            message(f'stale input {path}')
            missing = True
            continue

        mtime = os.path.getmtime(path)
        if newest_mtime == None or mtime > newest_mtime:
            newest_mtime = mtime

    if missing: return "missing"

    assert newest_mtime != None, 'invalid newest_mtime!'
    return newest_mtime

def check_output_mtimes(output):
    '''
    return oldest mtime if all outputs present
    return "missing" if any missing or stale
    return "empty" if job does not specify any outputs
    '''

    if output.is_empty(): return "empty"

    oldest_mtime = None
    missing = False

    for item,path in output.items():
        if not os.path.exists(path):
            message(f'missing output {path}')
            missing = True
            continue

        if is_stale(path):
            message(f'stale output {path}')
            missing = True
            continue

        mtime = os.path.getmtime(path)

        if oldest_mtime == None or mtime < oldest_mtime:
            oldest_mtime = mtime

    if missing: return "missing"

    assert oldest_mtime != None, 'invalid oldest_mtime!'
    return oldest_mtime

def generate_shell_commands(action,job_list,shell):

    for job_numb,job in enumerate(job_list):
        #check all inputs
        #returns newest mtime, "missing" or "empty"
        newest_input = check_input_mtimes(job['input'])

        #returns oldest mtime, "missing" or "empty"
        oldest_output = check_output_mtimes(job['output'])

        run = False

        #one or more inputs missing, job not runnable
        if newest_input == "missing":
            #flag job for removal from the list
            message(f'one or more inputs missing or stale, job not runnable')
            run = False

        #one or more outputs missing, job needs to run
        elif oldest_output == "missing":
            message(f'one or more outputs missing or stale, job needs to run')
            run = True

        elif 'run' in action and action['run'] == 'always':
            warning('"run:always" forcing job to run')
            run = True

        #no way to check output mtimes wrt input
        elif newest_input == "empty":
            if oldest_output == "empty":
                message('running due to no inputs or outputs being specified, use run:"never" to prevent')
                run = True

            #all outputs present, assume they are good
            else:
                message("no inputs specified, assuming existing outputs are good")
                run = False

        else: #newest_input has an mtime
            #no outputs specified and run:never already prevented reaching this point
            #in the code therefore run the job
            if oldest_output == "empty":
                message('running due to no outputs being specified, use run:"never" to prevent')
                run = True

            else: #oldest_output has an mtime
                if oldest_output < newest_input:
                    #outputs older than inputs
                    message('outputs older than inputs, running')
                    run = True
                else:
                    #all outputs look equal or newer than inputs
                    #we allow equal to count as newer
                    #so that a symlink to a file is not stale with respect to
                    #the target file so we can make a symlink to a file
                    #and not have it always look stale
                    #because the default "getmtime" of a symlink
                    #is that of the target not the actual symlink
                    message('outputs not older than inputs, assuming good')
                    run = False

        if run == False:
            #flag as defunct and move onto next potential job in the list
            job_list[job_numb] = None
            continue

        #job is going to run

        #remove stale output files/symlinks/directories before action is run
        handle_stale_outputs(action,job['output'])

        #create any missing output *parent* directories
        if action['ym/missing_parent_dir'] == 'create':
            create_output_dirs(action,job['output'])

    #remove non-runnable jobs
    job_list = [job for job in job_list if job is not None]
    shell_list = []

    for job_numb,job in enumerate(job_list):
        #merge input and output filenames and placeholders into action variables
        config = Conf(action)
        config.update(job['input'])
        config.update(job['output'])

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

    if not is_active(): return

    dt_old = datetime.datetime.fromisoformat('1980-01-01').timestamp()
    os.utime(path, (dt_old, dt_old))

def recycle_item(config,path):
    'move the path into the recycle bin'
    
    if not is_active(): return

    recycle_bin = config['ym/recycle_bin']

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
        if parent == '': continue #current working directory
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

def warning(item, end='\n', timestamp=True):
    message(col["yellow"] + warning_prefix + item + col["none"], end=end, timestamp=timestamp)

def error(item, end='\n', timestamp=True):
    message(col["red"] + error_prefix + item + col["none"], end=end, timestamp=timestamp)

def header(item,end='\n',timestamp=True):
    line = '-'*(70-len(item))
    item = f'{item}  {line}'
    blank()
    message(col["green"]+item+col["none"],end=end,timestamp=timestamp)
    flush()

def flush():
    message('',end='',timestamp=False)

def blank():
    message('',end='\n',timestamp=False)

def message(item,end='\n',timestamp=True):
    if timestamp:
        item = timestamp_now_nice() + ' ' + str(item)

    if meta['args'].quiet != True:
        display_message(item+end)
        sys.stdout.flush()

    if meta['args'].no_logs != True:
        path = os.path.join(meta['log_dir'],
                            meta['prefix']+meta['start_time']+'.messages')

        try:
            with open(path,'a') as f:
                f.write(item+end)
                f.close()
        except:
            print(timestamp_now_nice() + ' ' + warning_prefix+f'log file {path} not writeable')

def display_message(item):
    if not 'prev_message' in meta:
        meta['prev_message'] = None
        meta['prev_count'] = 0

    if item == meta['prev_message']:
        meta['prev_count'] += 1
        return

    if meta['prev_message'] != None:
        if meta['prev_count'] > 1:
            if meta['prev_message'].endswith('\n'):
                meta['prev_message'] = meta['prev_message'][:-1]
                meta['prev_message'] += f' [ x{meta["prev_count"]} ]\n'
            else:
                meta['prev_message'] += f' [ x{meta["prev_count"]} ]'

        print(meta['prev_message'],end='')

    meta['prev_message'] = item
    meta['prev_count'] = 1

def handle_stale_outputs(config,outputs):
    '''
    the job is going to be rerun therefore any outputs already present
    are treated as stale and recycled or deleted
    '''

    for item in outputs.keys():
        path = outputs[item]
        if not os.path.exists(path): continue

        if os.path.islink(path) or os.path.isfile(path):
            if config['ym/stale_output_file'] == 'delete':
                message(f'deleting stale output file {path}')
                remove_item(path)
            elif config['ym/stale_output_file'] == 'recycle':
                message(f'recycling stale output file {path}')
                recycle_item(config,path)
            elif config['ym/stale_output_file'] == 'ignore':
                message(f'ignoring stale output file {path}')
            elif config['ym/stale_output_file'] == 'stale':
                message(f'explicitly marking as stale output file {path}')
                make_stale(path)
            else:
                raise Exception(f'unknown option for stale_output_file: {config["ym/stale_output_file"]}')

        elif os.path.isdir(path):
            if config['ym/stale_output_dir'] == 'delete':
                message(f'deleting stale output directory {path}')
                remove_tree(path)
            elif config['ym/stale_output_dir'] == 'recycle':
                message(f'recycling stale output directory {path}')
                recycle_item(config,path)
            elif config['ym/stale_output_dir'] == 'ignore':
                message(f'ignoring stale output directory {path}')
            elif config['ym/stale_output_dir'] == 'stale':
                message(f'explicitly marking as stale output dir {path}')
                make_stale(path)
            else:
                raise Exception(f'unknown option for stale_output_dir: {config["ym/stale_output_dir"]}')
                
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
            if config['ym/failed_output_file'] == 'delete':
                message(f'deleting failed output file {path}')
                remove_item(path)
            elif config['ym/failed_output_file'] == 'recycle':
                message(f'recycling failed output file {path}')
                recycle_item(config,path)
            elif config['ym/failed_output_file'] == 'stale':
                message(f'marking as stale failed output file {path}')
                make_stale(path)
            elif config['ym/failed_output_file'] == 'ignore':
                message(f'ignoring failed output file {path}')
            else:
                raise Exception(f'unknown option for failed_output_file: {config["ym/failed_output_file"]}')

        elif os.path.isdir(path):
            if config['ym/failed_output_dir'] == 'delete':
                message(f'deleting failed output directory {path}')
                remove_tree(path)
            elif config['ym/failed_output_dir'] == 'recycle':
                message(f'recycling failed output directory {path}')
                recycle_item(config,path)
            elif config['ym/failed_output_dir'] == 'stale':
                message(f'marking as stale failed output directory {path}')
                make_stale(path)
            elif config['ym/failed_output_dir'] == 'ignore':
                message(f'ignoring failed output directory {path}')
            else:
                raise Exception(f'unknown option for failed_output_dir: {config["ym/failed_output_dir"]}')
                
        else:
            raise Exception(f'unsupported output type {path}')

def rmnl(item):
    '''
    remove one trailing newline if present
    because it will be added back when they get joined
    '''

    assert type(item) == str
    if item.endswith('\n'): return item[:-1]
    return item

def generate_full_command(config,shell):
    'generate the bash commands to run before the job commands'

    cmd_list = []

    if 'ym/bash_setup' in config and config['ym/bash_setup'] != '':
        cmd_list.append(rmnl(config['ym/bash_setup']))

    if 'conda' in config and config['conda'] != '':
        if 'ym/conda_setup' in config and config['ym/conda_setup'] != '':
            cmd_list.append(rmnl(config['ym/conda_setup']))
        cmd_list.append(rmnl(f'conda activate {config["ym/conda_prefix"]}{config["conda"]}'))

    cmd_list.append(rmnl(shell))

    return '\n'.join(cmd_list)

def add_env(config,env):
    'add key:value pairs from config["env"] to env'

    #no environment variables provided
    if not 'env' in config: return

    for key,value in config['env'].items():
        assert type(value) == str, f'non-string value found in env field {key}, all env fields must be simple strings'
        env[key] = value

def generate_job_environment(config,job_numb,njobs):
    'copy local environment with a few adjustments'

    env = copy.deepcopy(os.environ)

    #add in items from the special 'env' field if present
    add_env(config,env)

    env[ config['ym/job_count'] ] = f'{njobs}'
    env[ config['ym/job_number'] ] = f'{job_numb+1}' # convert from 0 to 1 based to match SGE_TASK_ID

    return env

def write_jobfile(action,shell_list):
    'save action config and shell_list as json'

    fnamebase = f'{action["ym/prefix"]}{timestamp_now()}.{action["name"]}'
    jobfile = os.path.join(meta['log_dir'],fnamebase+'.jobs')
    payload = {'action':action.getdict(),'shell_list':shell_list}

    with open(jobfile,'w') as f:
        json.dump(payload,f,sort_keys=False,indent=2)
        f.write('\n')

    return fnamebase

def colorize_command(cmd,c):
    '''
    put colour escape sequences at the start of each line
    to make less -SR work the same as direct output
    '''

    tmp = cmd.replace("\n",f"\n{c}")
    tmp = f"{c}{tmp}{col['none']}"

    return tmp

def execute_command(config,job_numb,cmd,env):
    'execute command locally'
    fname = f'{config["ym/prefix"]}{timestamp_now()}.{config["name"]}'
    foutname = os.path.join(meta['log_dir'],fname+'.out')
    ferrname = os.path.join(meta['log_dir'],fname+'.err')

    message(f'job {job_numb+1} executing locally...')

    if cmd.endswith('\n'): end = ''
    else:                  end = '\n'

    tmpcmd = colorize_command(cmd,col['green'])
    message(f'{tmpcmd}',end=end,timestamp=False)

    if not is_active():
        #dry-run: signal job completed ok without running it
        message(f'pipeline inactive, skipping actual command execution')
        return False

    flush()

    failed = False
    fout = open(foutname,'w')
    ferr = open(ferrname,'w')

    try:
        subprocess.run(cmd,env=env,shell=True,check=True,stdout=fout,stderr=ferr)
    except subprocess.CalledProcessError:
        failed = True

    fout.close()
    ferr.close()

    if failed: error('command failed')
    else:      message('command completed normally')

    flush()
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
    env["yamlmake"] = sys.argv[0]
    env["log_dir"] = meta["log_dir"]

    for line in f_in: f_out.write(line.format(**env))

    f_out.close()
    f_in.close()

def submit_job_qsub(action,shell_list,job_list):
    '''
    issue the qsub command to spawn an array of jobs
    wait for completion before checking the status of each job
    '''

    jobname = write_jobfile(action,shell_list)
    qsub_script = os.path.join(meta['log_dir'],jobname+'.qsub')
    jobfile = os.path.join(meta['log_dir'],jobname+'.jobs')
    njobs = len(shell_list)

    write_qsub_file(action,qsub_script,jobname,njobs,jobfile)

    cmd = f"qsub -sync y {qsub_script}"

    env = copy.deepcopy(os.environ)

    message(f'executing {cmd} with {njobs} task(s)...')

    something_failed = False

    if is_active():
        foutname = os.path.join(meta['log_dir'],jobname+'.out')
        ferrname = os.path.join(meta['log_dir'],jobname+'.err')
        fout = open(foutname,'w')
        ferr = open(ferrname,'w')
        #sys.stdout.flush()
        flush()

        try:
            subprocess.run(cmd,env=env,shell=True,check=True,stdout=fout,stderr=ferr)
        except subprocess.CalledProcessError:
            something_failed = True

        fout.close()
        ferr.close()

        if something_failed: warning(f'qsub call failed')
        else:                message(f'qsub call succeeded')

    #delay to allow for shared filesystem latency on status files
    message(f'sleeping for {action["ym/remote_delay_secs"]} seconds to allow for filesystem latency...')
    flush()
    if activity_state() == 'active': time.sleep(int(action['ym/remote_delay_secs']))

    #check each individual job's status file and expected outputs
    #some jobs may not have failed even if qsub returned an error code?
    for job_numb,item in enumerate(shell_list):
        status_file = f'{meta["log_dir"]}/{jobname}.{job_numb+1}.status'

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

        if not failed:
            #job has also failed if not all required outputs were created
            failed = verify_expected_outputs(action,job_list[job_numb]["output"],job_list[job_numb]["input"])

        if failed:
            #deal with output of failed job
            handle_failed_outputs(action,job_list[job_numb]["output"])
            something_failed = True

    flush()
    return something_failed

def qsub_execute_job(jobfile):
    '''
    runs on a compute node under grid engine control
    execute one job of a job array spawned by qsub
    invoked automatically by yamlmake using the --qsub option
    '''

    action,shell_list = read_jobfile(jobfile)
    action = Conf(action)

    assert jobfile.endswith('.jobs')

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
    if any are missing or stale flag job has failed

    implies jobs must use touch to update any output file that does not need altering
    or else exclude it from the list of required outputs

    return True to indicate failure
    '''

    if outputs.is_empty():
        message('no outputs specified')
        return False

    #find newest input file
    #return newest_mtime, "empty" or "missing"
    newest_mtime = check_input_mtimes(inputs)

    failed = False
    show = True
    for item in outputs.keys():
        path = outputs[item]
        if not os.path.exists(path):
            warning(f'missing output {path}')
            failed = True 

        elif is_stale(path):
            warning(f'stale output {path}')
            failed = True

        elif newest_mtime == 'empty':
            if show:
                message(f'cannot verify output freshness of {path} due to no inputs being specified')
                show = False

        elif newest_mtime == 'missing':
            #this will happen if the job deletes or moves any of its inputs
            if show:
                warning(f'cannot verify output freshness of {path} due to missing input(s)')
                show = False

        elif os.path.getmtime(path) < newest_mtime:
            #expected output file exists but appears stale
            #ie we assume it was not updated by the command
            message(f'stale output {path}')
            failed = True

    if failed:
        #job has failed
        return True

    #signal job seems to have completed okay
    if newest_mtime not in ('empty','missing'):
        message('all outputs seem fresh')

    return False

def execute_jobs(action,shell_list,job_list):
    something_failed = False

    if not 'exec' in action or action['exec'] == 'local':
        #local serial execution
        for job_numb,item in enumerate(shell_list):
            cmd = generate_full_command(action,item)
            env = generate_job_environment(action,job_numb,len(shell_list))

            #reports only explicit job failure from exit code
            failed = execute_command(action,job_numb,cmd,env)

            if not failed:
                #job has also failed if not all required outputs were created/updated
                failed = verify_expected_outputs(action,job_list[job_numb]["output"],job_list[job_numb]["input"])

            if failed:
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
        header(f'[{action["name"]}] skipped: pipeline inactive')
        return False

    if 'run' in action and action['run'] == 'never':
        header(f'[{action["name"]}] skipped: never run is set')
        return False

    if activity_state() == 'dryrun':
        header(f'[{action["name"]}] dry-running')
    else:
        header(f'[{action["name"]}] running')

    #generate list of potential jobs by filling out glob and list placeholders
    #this step does not pay any attention to file time stamps
    #but will reject anything with missing input file(s) already
    job_list = generate_job_list(action,input,output)
    message(f'{len(job_list)} potential job(s) after placeholder expansion')
    if len(job_list) == 0: return False

    #determine if all inputs are present and non-stale
    #determine if any outputs need (re)generating
    job_list,shell_list = generate_shell_commands(action,job_list,shell)
    if len(shell_list) == 0:
        message('no runnable job(s) after input/output file checking')
        return False

    message(f'==> {col["cyan"]}{len(shell_list)} runnable job(s){col["none"]} <==')

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
    process the YAML pipeline that was already loaded from path
    '''

    #initially set config to default values if none provided
    if config == None:
        #set conf to defaults
        config = Conf(src="defaults")

    #on first call store args and init stateful variables that
    #implement run-only, run-from and run-until options
    init_meta(args,config)

    counter = 0

    header(f"starting pipeline [{path}]")

    while counter < len(pipeline):
        item = pipeline[counter]
        counter += 1
        assert type(item) == dict
        assert len(item) == 1
        item_type = list(item.keys())[0]

        if item_type == 'action':
            if process_action(config,item[item_type],path):
                if activity_state() == 'dryrun':
                    warning('action failed due to dryrun mode, continuing pipeline anyway')
                else:
                    warning('aborting pipeline due to failed action')
                    return 2

        elif item_type == 'config':
            #add new config to the existing one, overriding any shared keys
            message(f'updating config')
            config.update(item[item_type])
            config.includes_and_loads(path)

        elif item_type == 'include':
            #toplevel include: load and insert the yaml items in place of the include item
            header(f'[{item[item_type]}] included')
            new_pipeline,new_path = load_pipeline(item[item_type],path) #new_path ignored
            counter -= 1
            del pipeline[counter]
            pipeline = pipeline[:counter] + new_pipeline + pipeline[counter:]

        elif item_type == 'module':
            #process a nested pipeline without affecting the config of any
            #following items
            header(f'[{item[item_type]}] module')
            new_pipeline,new_path = load_pipeline(item[item_type],path)
            sub_config = Conf(config)
            process(new_pipeline,new_path,config=sub_config)

        else:
            raise Exception(f"unsupported item type {item_type}")

    message(f"reached end of pipeline {path}")
    flush()

    return 0
