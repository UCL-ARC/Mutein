'''
------------------------
RSA 29/03/22
------------------------
Helper module for pipeline script for foldx job on Myriad
----
'''

import os

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
    return params

def configparams(filename):    
    #set up some defaults for any batch to run without paramaters
    params = {}
    params['jobs'] = '1234567'
    params['pdb'] = '6vxx'
    params['name'] = '6vxx_50'
    params['split'] = '50'
    params['mutation'] = '.'
    params['time'] = '.'
    params['combos'] = '63'
    params['variant'] = 'Alpha'    
    if filename != '':
        with open(filename) as fr:
            cfgcontent = fr.readlines()
            for line in cfgcontent:
                line = line.strip()
                params = addlinetoparams(line,params)
    if 'variantfile' not in params:
        params['variantfile'] = params['pdb'] + '_vars'

    #Code to decide if it is test or live environment can be overridden from the command line
    dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'        
    if 'wsl$' in dir_path or '/rachel/' in dir_path:# Possibly need to change to rlevant python exe and path for local runs
        env = 'python'
    else:
        env = 'hpc'
    params['env'] = env

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







