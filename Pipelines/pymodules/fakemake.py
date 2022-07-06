import yaml
import re
import copy
import json
import os
import glob

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
#= file glob: creates separate jobs
###* file glob: creates file list within single job <= not implemented
reserved_chrs = ['$','=']       #,'*']

#default conda and exec values
default_conda = 'fm_default_env'
default_exec = 'local'

# def _sub_vars(config,extra=None):
#     'substitute all simple placeholders or fail trying'
#     counter = 10
#     while True:
#         counter -= 1
#         finished = True
#         for key,value in config.items():
#             for m in re.finditer(ph_regx,value):
#                 name = m.group(0)[1:-1]
#                 #if name[0] in reserved_chrs: continue #ignore special placeholders

#                 if extra == None or name in config:
#                     new_value = config[name]
#                 else:
#                     new_value = extra[name]
#                 config[key] = value[:m.start(0)] + new_value + value[m.end(0):]
#                 finished = False
#         if finished: break
#         assert counter > 0, 'unable to resolve all placeholders'
        
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
    #show(job_list,"jobs")

    #determine if all inputs are present
    #and if outputs need to be regenerated
    validate_jobs(rule,job_list,shell)

def process(config,pipeline):
    #pipeline yaml must be a list of items all called "rule"
    for i,item in enumerate(pipeline):
        assert type(item) == dict
        assert len(item) == 1
        item_type = list(item.keys())[0]
        assert item_type == 'rule'
        #show(item[item_type],item_type)
        process_rule(config,item[item_type])

def setup_rule(config,rule):
    assert type(rule) == dict

    #insert any missing defaults
    if not 'exec' in rule: rule['exec'] = default_exec
    if not 'conda' in rule: rule['conda'] = default_conda

    #check all keys and values are simple strings
    #check for forbidden keys
    for key,value in rule.items():
        assert key not in nonrule_keys
        assert type(key) == str
        if key not in ['input','output']:
            #general rule options must be simple strings
            assert type(value) == str
        else:
            #input and output can also be dictionaries of strings
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

    #merge temporary copy of config into the rule, giving rule priority
    _config = copy.deepcopy(config)
    _config.update(rule)
    rule.update(_config)

    #substitute any environment variables
    sub_environ(rule)
    sub_environ(input)
    sub_environ(output)
    sub_environ(shell)

    #substitute any fakemake variables except for shell
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

    #convert simple form input/output into ordered list form
    if type(input) == str: input = {'input':input}
    if type(output) == str: output = {'output':output}
    shell = {'shell':shell}

    return input,output,shell

def setup_config(pipeline):
    #set any defaults here
    config = {'example_default':'example_value'}

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

def validate_jobs(rule,job_list,shell):
    for job in job_list:
        job['shell'] = copy.deepcopy(shell)

        #merge input and output filenames into rule variables
        config = copy.deepcopy(rule)
        config.update(job['input'])
        config.update(job['output'])

        #substitute remaining placeholders in shell command
        sub_vars(job['shell'],config)

        show(job['shell'],"shell")
