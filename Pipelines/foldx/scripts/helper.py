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
environments['empty_hpc'] = ['foldx', 'python','empty_hpc'] #just prints out what it would run
environments['empty_python'] = ['foldx', 'python','empty_python'] #just prints out what it would run
environments['python'] = ['foldx', 'python','python'] #just prints out what it would run
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
    if '=' in arg and "@" not in arg:
        args = arg.split('=')
        p,v = args[0],args[1]
        params[p]=v    
    return params

def addpipelinetoparams(arg,params):    
    args = []
    if '=' in arg and "@" in arg:
        args = arg.split('@')
        print(args)
        id = args[0].split('=')[1]  
        pv = args[1].split('=')
        p,v = str(pv[0]),str(pv[1])
        print(args[0],args[1],id,p,v)
        if id not in params:
            params[id] = {}
        params[id][p]=str(v)    
    return params

def configparams(pdb):    
    #set up some defaults for any batch to run without paramaters    
    params = {}
    params['jobs'] = '1234567'
    params['chain'] = 'A'
    params['pdb'] = '6vxx'
    params['name'] = '6vxx_50'    
    params['row'] = '1'
    params['mutation'] = '.'
    params['time'] = '.'    
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

def configpipelineparams(pdb):    
    #set up some defaults for any batch to run without paramaters    
    params = {}    
    if pdb != '':
        dir_path = os.path.dirname(os.path.realpath(__file__))
        dir_path = dir_path[:-7]
        input_path = dir_path + 'inputs/'         
        configfile = input_path + pdb + '/config.cfg'
        with open(configfile) as fr:
            cfgcontent = fr.readlines()
            for line in cfgcontent:
                line = line.strip()
                params = addpipelinetoparams(line,params)        
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
            
def pipelineparams(argvs,params):            
    for i in range(1,len(argvs)):
        arg = argvs[i]
        params = addpipelinetoparams(arg,params)                
    return params

def mergeparams(configparams, jobparams):
    #the job params take precendence
    for cfg,val in jobparams.items():
        configparams[cfg] = val
    return configparams

def add_remove_path_levels(path,levels=0,dir=''):
    dirs = path.split('/')
    if levels > 0:
        dirs = dirs[:-1*levels]
    if dir != '':
        dirs.append(dir)
    retpath = ''
    for d in dirs:
        retpath += d + '/'
    if not os.path.exists(retpath):
        os.mkdir(retpath)    
    return retpath
    
def get_make_paths(pdb,name):
    paths_dic = {}
    # The path structure is relative to.....it could be a given path
    # For now the path structure is relative to the script file
    # (But it will be easy to change)    
    dir_script_path = os.path.dirname(os.path.realpath(__file__))
    dir_path = add_remove_path_levels(dir_script_path,1)
    inputs_path = add_remove_path_levels(dir_path,0,'inputs')
    interim_path = add_remove_path_levels(dir_path,0,'interim')
    thruputs_path = add_remove_path_levels(dir_path,0,'thruputs')
    outputs_path = add_remove_path_levels(dir_path,0,'outputs')
    inputs_path = add_remove_path_levels(inputs_path,0,pdb)
    thruput_path = add_remove_path_levels(thruputs_path,0,pdb)
    interim_path = add_remove_path_levels(interim_path,0,name)
    output_path = add_remove_path_levels(outputs_path,0,name)            
    return inputs_path, thruput_path, interim_path, output_path

def goto_job_dir(dir_path,args,params,name):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)  
    os.chdir(dir_path)
    inputs_file = name+'.log'
    with open(inputs_file, 'w') as fw:
        for arg in args:
            fw.write(str(arg) + ' ')
        fw.write('\n')
        for cfg,val in params.items():
            fw.write(str(cfg) + '=' + str(val) + '\n')
            






