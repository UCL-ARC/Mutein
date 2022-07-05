import yaml
import re
import copy
import json
import os
import glob

def sub_place_holders(key,avail,config,pattern):
    for m in re.finditer(pattern,config[key]):
        placeholder = m.group(0)[1:-1]
        if placeholder.startswith('_'): continue

        if placeholder.startswith('$'):
            new_value = os.environ[placeholder[1:]]
        else:
            if not placeholder in avail: continue
            new_value = config[placeholder]

        config[key]  = config[key][:m.start(0)] + new_value + config[key][m.end(0):]
        
def pending(value,pattern):
    if re.search(pattern,value) != None: return True
    return False

def substitute(config):
    pattern = config['__pattern']
    avail = set()
    pend = set()

    for key,value in config.items():
        if key.startswith('_'):
            continue               #protected core setting or wildcard
        elif pending(value,pattern):
            pend.add(key)      #needs substituting still
        else:
            avail.add(key)     #final value alrady presenrt

    for key in pend:
        sub_place_holders(key,avail,config,pattern)

def show(item,label,indent=2):
    print(label)
    print(json.dumps(item,sort_keys=True,indent=indent))
    print()

def process_rule(rule,config):

    #apply place holder substitution
    substitute(config)

    #check option validity
    assert config['exec'] in ['local','qsub']

    show(config,"rule")

def process(pipeline):
    #load some defaults
    config = {'__pattern':'\{[^\}]*\}'}

    for i,item in enumerate(pipeline):
        assert type(item) == dict
        assert len(item) == 1
        key = list(item.keys())[0]

        if key == 'config':
            if i != 0:
                raise Exception("config must be first item in pipeline file")

            config.update(item[key])
            show(config,"config")

        elif key == 'rule':
            process_rule(item[key],config)

def parse_yaml(fname):
    'parse yaml into list of items'
    with open(fname) as f:
        result = yaml.safe_load(f)

    assert type(result) == list
    
    return result

#def find_env_variables(item):

# test_input = "datasets/{=dataset}/{=subset}/{=accession}/{accession}_{*read}.fastq.gz"
# test_output = "datasets/{dataset}/{subset}/{accession}/{accession}.info"
#in general it would be lists of input and output files not just one or each
#"datasets/(?P<dataset>[^/]*)/(?P<subset>[^/]*)/(?P<accession>[^/]*)/(?P=<accession>)"
def expand_io_paths(input,output):
    pattern = '\{[^\}]+\}'
    input_glob = ''
    input_regx = '^'
    pholders = []
    ditems = {}
    prev_end = 0

    for m in re.finditer(pattern,input):
        name = m.group(0)[1:-1]
        start = m.start(0)

        input_glob += input[prev_end:start] + '*'
        input_regx += input[prev_end:start]

        if name.startswith('='):
            #placeholder defines a set of separate *jobs*
            name = name[1:]
            assert name not in ditems
            pholders.append(name)
            ditems[name] = 'job'
            input_regx += '(?P<' + name + '>[^/]+)'
        elif name.startswith('*'):
            #placeholder defines a *list* of filenames within a job
            name = name[1:]
            assert name not in ditems
            pholders.append(name)
            ditems[name] = 'list'
            input_regx += '(?P<' + name + '>[^/]+)'
        else:
            #back reference to previous job or list placeholder
            assert name in ditems
            input_regx += '(?P=' + name + ')'

        prev_end = m.end(0)

    input_glob += input[prev_end:]
    input_regx += input[prev_end:] + '$'

    print(input)
    print(input_glob)
    print(input_regx)

    for path in glob.iglob(input_glob):
        print(path)
        m = re.fullmatch(input_regx,path)
        if m is not None:
            print(m.group(0))
            for name in pholders: print(m.group(name))
