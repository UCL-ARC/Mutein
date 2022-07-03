import yaml
import re
import copy
import json
import os

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
def expand_io_paths(input,output):
    pattern = '\{[^\}]*\}'
    input_glob = ''
    items = []
    for m in re.finditer(pattern,input):
        name = m.group(0)[1:-1]
        if name.startswith('='):
            ptype = 'job'
        elif name.startswith('*'):
            ptype = 'list'
        else:
            ptype = 'normal'

        if len(items) > 0:
            start = items[-1][3]
        else:
            start = 0
        end = m.start(0)
        input_glob += input[start:end]

        input_glob += '*'
        items.append([name,ptype,m.start(0),m.end(0)])

    start = items[-1][3]
    input_glob += input[start:]

    print(input)
    print(input_glob)
