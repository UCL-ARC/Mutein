'''
------------------------
RSA 29/03/22
------------------------
Helper module for pipeline script for foldx job on Myriad
----
'''

import os

##################################### Environments exes and paths  ##########################################################
######################################### !!!!! TODO !!!!! ##################################################################
## New users, add your environment details here to and check them in to avoid having to mess about locally all the time ##
# For example, when I run locally, I want to use python instead of hpc, and my exe path to foldx and python is different (environment variable incompetence)
environments = {}
environments['empty'] = ['foldx', 'python','empty'] #just prints out what it would run
environments['CI'] = ['~/UCL/libs/foldx5/foldx', '/bin/python3','python'] #continuous integration
environments['myriad'] = ['foldx', 'python','hpc']
environments['myriad_tst'] = ['foldx', 'python','python']
environments['rachel'] = ['~/UCL/libs/foldx5/foldx', '/bin/python3','python']

def getenvironment(user=''):
    '''Automatically recognise the environment though to can be overridden by exploicitly passing it in
    Returns: environment, tuple(foldx path, python path, hpc or python)
    '''    
    if user == '':
        dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'  
        if '/rachel/' in dir_path:
            user = 'rachel'
        #elif add your user env here to default
        else:
            user = 'myriad'
    if user in environments:        
        return user,environments[user]    

    return None,None

##################################### Environments exes and paths  ##########################################################

### These functions consistently handle the paramater inputs for the script, merging config and overrides
def addlinetoparams(arg,params):    
    args = []
    if '=' in arg:
        args = arg.split('=')
    if 'jobs=' in arg:
        params['jobs'] = args[1]
    elif 'pdb=' in arg:
        params['pdb'] = args[1]
    elif 'name=' in arg:
        params['name'] = args[1]
    elif 'split=' in arg:
        params['split'] = args[1]
    elif 'mutation=' in arg:
        params['mutation'] = args[1]
    elif 'time=' in arg:
        params['time'] = args[1]
    elif 'variant=' in arg:
        params['variant'] = args[1]
    elif 'variantfile=' in arg:
        params['variantfile'] = args[1]
    elif 'combos=' in arg:
        params['combos'] = args[1]
    elif 'env=' in arg:
        params['env'] = args[1]
    elif 'user=' in arg:
        params['user'] = args[1]
    elif 'chain=' in arg:
        params['chain'] = args[1]
    elif 'row=' in arg:
        params['row'] = args[1]
    return params

def configparams(pdb):    
    #set up some defaults for any batch to run without paramaters    
    params = {}
    params['jobs'] = '1234567'
    params['chain'] = 'A'
    params['pdb'] = '6vxx'
    params['name'] = '6vxx_50'
    params['split'] = '50'
    params['row'] = '1'
    params['mutation'] = '.'
    params['time'] = '.'
    params['combos'] = '63'
    params['variant'] = 'Alpha'        
    if pdb != '':
        dir_path = os.path.dirname(os.path.realpath(__file__))
        dir_path = dir_path[:-7]
        input_path = dir_path + 'inputs/'         
        configfile = input_path + pdb + '/config.cfg'
        with open(configfile) as fr:
            cfgcontent = fr.readlines()
            for line in cfgcontent:
                line = line.strip()
                params = addlinetoparams(line,params)
    if 'variantfile' not in params:
        params['variantfile'] = params['pdb'] + '_vars'

    #Code to decide if it is test or live environment can be overridden from the command line
    user, envs = getenvironment()    
    params['user'] = user
    return params

def inputparams(argvs):        
    params = {}
    for i in range(1,len(argvs)):
        arg = argvs[i]
        params = addlinetoparams(arg,params)            
    if 'pdb' not in params:
        params['pdb'] = '6vxx'    
    if 'configfile' not in params:
        params['configfile'] = '../inputs/'+ params['pdb'] + '/config.cfg'    
    return params
            

def mergeparams(configparams, jobparams):
    #the job params take precendence
    for cfg,val in jobparams.items():
        configparams[cfg] = val
    return configparams

def get_make_paths(pdb,name):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    dir_path = dir_path[:-7]
    if not os.path.exists(dir_path + 'interim/'):
        os.mkdir(dir_path + 'interim/')
    if not os.path.exists(dir_path + 'thruputs/'):
        os.mkdir(dir_path + 'thruputs/')
    if not os.path.exists(dir_path + 'outputs/'):
        os.mkdir(dir_path + 'outputs/')
    input_path = dir_path + 'inputs/' + pdb + '/'
    thruput_path = dir_path + 'thruputs/' + pdb + '/'
    interim_path =     dir_path + 'interim/' + name + '/'
    output_path =     dir_path + 'outputs/' + name + '/'
    if not os.path.exists(thruput_path):
        os.mkdir(thruput_path)
    if not os.path.exists(interim_path):
        os.mkdir(interim_path)
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    return input_path, thruput_path, interim_path, output_path



