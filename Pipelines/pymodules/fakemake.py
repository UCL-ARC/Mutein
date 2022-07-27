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

#regex patterns to match non nested {placeholders} and {$environment variables}
environ_regx  = r'\{\$.+?\}' #{$name}
scalar_regx   = r'\{%.+?\}'  #{%scalar}
list_regx     = r'\{&.+?\}'  #{&list[]}
glob_regx     = r'\{\*.+?\}' #{*glob} ==> job scalar
listjob_regx  = r'\{=.+?\}'  #{=list} ==> job scalar

#required and forbidden keys (after input/output/shell are removed from rule)
required_keys = {
    "config": [  ],
    "rule":   [ "name" ],
    "input":  [  ],
    "output": [  ],
    "shell":  [ "shell" ],
}

forbidden_keys = {
    "config": [ 'input','output','name','rule','config','module','include','import' ],
    "rule":   [ "config" ],
    "input":  [  ],
    "output": [  ],
    "shell":  [  ],
}

allow_lists = {
    "config": True,
    "rule":   True,
    "input":  False,
    "output": False,
    "shell":  False,
}

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

class Conf:
    def __init__(self,src=None):
        #note both scalars and lists should be ordered dictionaries
        #(requires python 3.7 or later)

        #config values that are all simple strings
        self.scalars = {}

        #config values that are lists of strings
        self.lists = {}

        if src == None:
            #initialise to empty
            pass
        elif src == "defaults":
            #initialise to defaults
            self.update(default_global_config)

        elif type(src) == Conf or type(src) == dict:
            #copy another config or dict
            self.update(src)

        else:
            raise Exception(f"unsupported source class {type(src)}")

    def __getitem__(self,key):
        'default getitem looks only in scalars'
        return self.scalars[key]

    def __setitem__(self,key,value):
        'default setitem only sets scalars'
        self.scalars[key] = value

    def show(self,label,indent=2):
        print(label)
        print(json.dumps(self.scalars,sort_keys=False,indent=indent))
        print(json.dumps(self.lists,sort_keys=False,indent=indent))
        print()

    def items(self):
        return self.scalars.items()

    def keys(self):
        return self.scalars.keys()

    def getlist(self,key):
        return self.lists[key]

    def setlist(self,key,value):
        self.lists[key] = value

        #self.check()

    def listitems(self):
        return self.lists.items()

    def update_dict(self,src_dict):
        '''
        add all key:value pairs from a raw yaml dictionary
        values of src overwrite any shared keys in self
        check for validity and split into scalars and lists
        '''

        assert type(src_dict) == dict

        for key,value in src_dict.items():
            assert type(key) == str

            if type(value) == str:
                self.scalars[key] = value
            elif type(value) == list:
                for x in value: assert type(x) == str
                self.lists[key] = value.copy()
            else:
                raise Exception(f'unsupported config item type {type(value)}')

    def update_conf(self,src):
        'add all key:value pairs from src, overwriting any shared keys'

        assert type(src) == Conf

        self.scalars.update(src.scalars)

        for key,value in src.lists.items():
            self.lists[key] = value.copy()

    def update(self,src):
        '''
        add all key:value pairs from another Conf object or dict
        values of src take priority (overwrite) any duplicate keys in self
        '''

        if type(src) == Conf:
            self.update_conf(src)
        elif type(src) == dict:
            self.update_dict(src)
        else:
            raise Exception(f"unsuported object type {type(src)}")

    def override(self,src_config):
        '''
        add all key:value pairs from another conf object
        values of self take priority (are retained) for any duplicate keys
        '''

        assert type(src_config) == Conf

        tmp_scalars = copy.deepcopy(src_config.scalars)
        tmp_scalars.update(self.scalars)
        self.scalars = tmp_scalars

        for key,value in src_config.lists.items():
            if not key in self.lists:
                self.lists[key] = value.copy()

    def sub_input_output(self,src):
        '''
        substitute only {=...} and {*...} placeholders
        or throw exception
        '''

        self.sub_values(src,glob_regx)
        self.sub_values(src,listjob_regx)

    def sub_vars2(self,src=None):
        '''
        substitute all placeholders except {=...} and {*...} 
        or throw exception
        '''

        #where to get the replacement values from
        if src == None: src = self

        counter = 10
        while True:
            changed = False
            counter -= 1
            if self.sub_values(os.environ,environ_regx): changed = True
            if self.sub_values(src,scalar_regx):         changed = True
            if self.sub_values(src,list_regx):           changed = True
            if not changed: break
            assert counter > 0, 'unable to resolve all placeholders in 10 iterations'

    def sub_values(self,src,regx):
        '''
        substitute placeholders
        ie {$environ}, {%scalars} or {&lists[i]}
        '''

        changed = False
        for key,value in self.scalars.items():
            self.scalars[key],flag = self.do_sub2(src,regx,value)
            if flag: changed = True

        for key,vlist in self.lists.items():
            for i,value in enumerate(vlist):
                vlist[i],flag = self.do_sub2(src,regx,value)
                if flag: changed = True

        return changed

    def do_sub2(self,src,regx,value):
        '''
        sub all simple placeholders in value from src
        src can be a dict (os.environ) or a Conf (ie Conf.scalars via __getitem__)
        '''

        changed = False

        while True:
            m = re.search(regx,value)
            if m == None: break #no more matches

            if m.group(0)[1] != '&':
                #not a list type placeholder
                name = m.group(0)[2:-1]   #{%name}, {$name}, {*name} or {=name}
                sub = src[name]
            else:
                #list type placeholder
                tmp_name = m.group(0)[2:-1]  # {&name[*]}
                assert tmp_name[-1] == ']'
                ind = tmp_name.index('[')
                name = tmp_name[:ind]        # name[*]
                subscript = tmp_name[ind+1:-1]

                if subscript in ['N']:
                    sub = str(len(src.getlist(name)))
                elif self.validate_subscript(subscript,len(src.getlist(name))):
                    sub = src.getlist(name)[int(subscript)]
                elif len(subscript) == 1:
                    if subscript == 't':   sub = '\t'.join(src.getlist(name))
                    elif subscript == 'n': sub = '\n'.join(src.getlist(name))
                    else:                  sub = subscript.join(src.getlist(name))
                    print(sub)
                else:
                    raise Exception(f"invalid list placeholder {tmp_name}")

            value = value[:m.start(0)] + sub + value[m.end(0):]
            changed = True

        return value,changed

    def validate_subscript(self,subscript,list_length):
        try:
            i = int(subscript)
        except:
            return False

        if i < -list_length or i >= list_length: return False

        return True

    # def sub_globs(self,extra):
    #     'substitute in any globbing placeholders from extra'
    #     return self.sub_pholders(extra,'=',gl_regx)

    def make_log_dir(self):
        if not os.path.exists(self['log_dir']):
            os.makedirs(self['log_dir'])

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
        #set conf to defaults
        config = Conf(src="defaults")
        config.show("default config")

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
            config.update(item[item_type])
            config.show("loaded config")

        elif item_type == 'load':
            #read in a list variable from a text file
            load_list(config,item[item_type])
            config.show("loaded list")

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
            sub_config = Conf(config)
            process(new_pipeline,new_path,config=sub_config)

        #show(config,"config")

def load_list(config,item):
    'load a list of values from a text file'
    
    name: "score"
    file: "filename"
    comment: "#"
    header: "0"
    sep: ","
    row: "2" | col: "score"

def load_pipeline(path,parent_file):
    'path is relative to the parent path'

    parent_path = os.path.dirname(parent_file)
    full_path = os.path.join(parent_path,path)

    assert full_path != parent_file
    return parse_yaml(full_path),full_path

def find_duplicates(conf_list):
    all_keys = set()

    for config in conf_list:
        for key in config.keys():
            assert not key in all_keys, f"duplicate key {key}"
            all_keys.add(key)

def setup_rule(config,rule):
    #separate out and canonicalise input, output and shell
    rule,input,output,shell = split_rule(rule)

    #merge config into rule but rule has priority
    rule.override(config)

    #try to ensure placeholder subtitution is not ambiguous
    find_duplicates([rule,input,output])

    #substitute all placeholders except for shell
    #which likely contains input/output variables yet to be determined
    rule.sub_vars2()
    input.sub_vars2(src=rule)
    output.sub_vars2(src=rule)

    #create log dir if missing
    rule.make_log_dir()

    return rule,input,output,shell

def split_rule(rule):
    '''
    separate out input, output and shell from rule
    convert all to Conf
    '''

    input = rule['input']
    output = rule['output']
    shell = rule['shell']

    del rule['input']
    del rule['output']
    del rule['shell']

    #convert simple form input/output into dictionary form
    if type(input) == str: input = {'input':input}
    if type(output) == str: output = {'output':output}
    shell = { 'shell':shell }

    rule = Conf(rule)
    input = Conf(input)
    output = Conf(output)
    shell = Conf(shell)

    return rule,input,output,shell

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
    glob or list expand the first input pattern to final file names if required
    fill out remaining input and output patterns with the matched names
    '''

    #the first input pattern triggers to glob/list expansion
    # primary_key = list(input.keys())[0]
    # primary = input[primary_key]

    #determine whether there are glob {*name} or list {=name} placeholders
    # mlist = re.search(listjob_regx,primary)
    # mglob = re.search(glob_regx,primary)

    # #cannot have both in same input pattern
    # assert not (mlist and mglob), "cannot use both globbing and list expansion together"

    #seed job "list" containing just the unexpanded input(s) and output(s)
    #which may have placeholders in need of expansion
    job_list = [{"input":input,"output":output,"matches":{}}]

    #expand input pattern lists
    new_list = []
    for job in job_list: new_list += generate_list_jobs(rule,job)
    job_list = new_list

    if len(job_list) == 0: return []

    print("after generate_list_jobs")
    for job in job_list:
        job['input'].show('job input')
        job['output'].show('job output')
        job['matches'].show('job matches')

    #expanded input pattern globs
    new_list = []
    for job in job_list: new_list += generate_glob_jobs(rule,job)
    job_list = new_list

    if len(job_list) == 0: return []

    print("after generate_glob_jobs")
    for job in job_list:
        job['input'].show('job input')
        job['output'].show('job output')
        job['matches'].show('job matches')

    exit()

    return job_list

def generate_list_jobs(rule,job):
    'expand job list to one item per input file pattern match (combo)'

    #identify all listjob type placeholders {=name} in input patterns
    ph_names = {}
    for key,value in job["input"].items():
        for m in re.finditer(listjob_regx,value):
            name = m.group(0)[2:-1]   #"{=name}" ==> "name"
            if name not in ph_names: ph_names[name] = True

    #nothing to expand if no list patterns present in input patterns
    if len(ph_names) == 0: return [ job ]

    #expand input and output patterns into complete paths
    #(excepting any remaining glob placeholders)
    #filled out with the values from the list(s)
    #where multiple lists are present expand to all combinations of the lists
    new_list = []
    meta_list = [rule.getlist(key) for key in ph_names]

    #product implements nested for loops over all the lists
    for value_list in itertools.product(*meta_list):
        src = { name:value_list[i] for i,name in enumerate(ph_names.keys()) }

        new_input = Conf(job["input"])
        new_output = Conf(job["output"])

        #substitute in the correct values from the list variables
        new_input.sub_values(src,listjob_regx)
        new_output.sub_values(src,listjob_regx)

        #add new matches to any existing ones
        new_matches = Conf(job["matches"])
        new_matches.update(src)
        new_list.append({"input":new_input,"output":new_output,"matches":new_matches})

    return new_list

def build_glob_and_regx(pattern):
    input_glob = ''
    input_regx = '^'
    ph_dict = {}
    prev_end = 0

    #convert glob-type placeholders into proper glob and regex query formats
    for m in re.finditer(glob_regx,pattern):
        name = m.group(0)[2:-1]   #{*name}
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

        prev_end = m.end(0)

    input_glob += pattern[prev_end:]
    input_regx += re.escape(pattern[prev_end:]) + '$'
    
    return input_glob,input_regx

def generate_glob_jobs(rule,job):
    'expand job list to one item per input file pattern match (combo)'

    glob_matches = {}

    #glob each input file pattern separately
    for key,pattern in job["input"].items():
        print(key,pattern)

        #convert the pattern-with-placeholders into a glob and regex string
        input_glob,input_regx = build_glob_and_regx(pattern)

        print(input_glob,input_regx)

        glob_matches[key] = []
        
        #find paths using iglob, match to placeholders using regex
        for path in glob.iglob(input_glob):
            m = re.fullmatch(input_regx,path)
            assert m is not None,"regex cannot extract placeholders from the globbed path!"

            print(path,m.group(0))

            #store the match as the completed path plus the placeholder values
            glob_matches[key].append([key,path,m.groupdict()])
        
        #no jobs are generated if any input pattern has zero matches
        if len(glob_matches[key]) == 0:
            print("no matches")
            return []

    #find glob_matches where the placeholders
    #have matching values across all input patterns
    #reject all others
    new_list = []
    meta_list = [ glob_matches[key] for key in glob_matches ]

    #product implements nested for loops over all the lists in meta_list
    for value_list in itertools.product(*meta_list):
        #build the potential new input and matches objects
        #check any shared placeholder values match
        #across all input patterns
        matches = {}
        inputs = {}

        conflict = False
        for key,path,groupdict in value_list:
            for name in groupdict:
                if name not in matches:
                    matches[name] = groupdict[name]
                elif matches[name] != groupdict[name]:
                    conflict = True
                    break
            if conflict: break
            inputs[key] = path

        #mismatching combo of input paths therefore no job generated
        if conflict: continue

        new_input = Conf(inputs)

        #substitute in the correct values from the matches
        new_output = Conf(job["output"])
        new_output.sub_values(matches,glob_regx)

        #add new matches to any existing ones from listjob generator
        new_matches = Conf(job["matches"])
        new_matches.update(matches)

        new_list.append({"input":new_input,"output":new_output,"matches":new_matches})

    return new_list

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
        #merge input and output filenames and placeholders into rule variables
        config = copy.deepcopy(rule)
        config.update(job['input'])
        config.update(job['output'])
        config.update(job['matches'])

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
    rule,input,output,shell = setup_rule(config,rule)
    rule.show("setup rule")

    print(f'\nrule: {rule["name"]}')

    #fill in any globbing placeholders in inputs and outputs
    #return list of values, one per potential job
    job_list = generate_job_list(rule,input,output)

    if len(job_list) == 0: return

    #determine if all inputs are present
    #and if outputs need to be regenerated
    job_list,shell_list = generate_shell_commands(rule,job_list,shell)

    print(f'{len(shell_list)} jobs generated')

    if len(shell_list) == 0: return

    #check current working directory agrees with configured value
    check_cwd(rule)

    #execute jobs locally or remotely from the jobfile
    execute_jobs(rule,shell_list,job_list)
