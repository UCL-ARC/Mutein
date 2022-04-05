'''
------------------------
RSA 29/03/22
------------------------
Pipeline overview script for foldx job on Myriad

This script is a parent script that runs the entire foldx pipeline at a top level.
It manages the dependencies beween the scripts and allows you to run as either HPC, python or empty mode

The scripts dependency is:

                                    SCRIPT 01: FOLDX REPAIR
                                    |                |
                        SCRIPT 02: SPLIT     SCRIPT 05: VARIANT SPLIT
                                |                       |
        [ARRAY JOBS]    SCRIPT 03: FOLDX POSSCAN     SCRIPT 06: FOLDX BUILD
                                |                       |
                        SCRIPT 04: AGGREGATE        SCRIPT 07: VARIANT AGGREGATE
------------------------
'''
import os
import subprocess
import helper as hlp
import pandas as pd

##### INPUTS #############################################
## Pipeline jobs sequence
# 1= repairing pdb
# 2= making param file for splits
# 3= performing the position scan ddg mutations (parallel)
# 4= aggregating 3
# 5= making params file for variants
# 6 = performing variant ddg (parallel)
# 7 = aggregating 6

def run_pipeline00(args):
    print('#### FOLDX PIPELINE - batch creation ####')
    ### Change into script directory
    dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
    print('## ... changing directory to',dir_path)
    os.chdir(dir_path)
    ### Process paramaters in order of preference, job, config, pipeline
    jobparams = hlp.inputparams(args)    
    pdb = ''
    if 'pdb' in jobparams:
        pdb = jobparams['pdb']
    cfgparams = hlp.configparams(pdb)    
    cfgplparams = hlp.configpipelineparams(pdb)
    print(cfgplparams)
    pipelineparams = hlp.pipelineparams(args,cfgplparams)
    print('Pipelines=',pipelineparams)
    runparams = hlp.mergeparams(cfgparams, jobparams) 
    #runparams['jobs'] = '1' ## if needed for testing purposes
    print('Params=',runparams)
    user = runparams['user']    
    user, (foldxe, pythonexe, env) = hlp.getenvironment(user)
    print('Environment=',user, foldxe, pythonexe,env)            
    # script extension is either sh or py depending on bash or python environment
    ext = '.sh'
    if env == 'python':
        ext = '.py'

    # The batch is defined in the file batch.csv    
    batch_df = pd.read_csv('batch.csv')
    batch_dic = {}
    for i in range(len(batch_df.index)):
        id = batch_df['id'][i]                
        script = batch_df['script'][i]                
        time = batch_df['time'][i]                
        dependency = batch_df['dependency'][i]                
        array = batch_df['array'][i]     
        batch_dic[str(id)] = (script,time,dependency,array)    
    #############################################################        
    dependencies = {}
    runs = []        
    for j in runparams['jobs']:        
        if str(j) in batch_dic:
            script,time,dependency,array = batch_dic[str(j)]
            #check there are no overrides from inuts
            if j in pipelineparams:
                idparams=pipelineparams[j]
                if 'time' in idparams:
                    time = idparams['time']
                if 'array' in idparams:
                    array = idparams['array']
            dep = "-1"
            if str(dependency) != "-1" and str(dependency) in runparams['jobs']:
                dep = dependency            
            runs.append([j,'qsub',script + ext,dep,time,array])            
            dependencies[str(j)] = str(dep)
                                                
    for job,exe,script,dependency,time,array in runs:
        #print(job,exe,script,dependency)
        if env == 'hpc':
            os.system('chmod +x ' + script)
        args = []
        if 'python' in env:
            args.append(pythonexe)        
        else:
            args.append('qsub')            
            if str(dependency) != "-1":
                dep = dependencies[int(dependency)]
                args.append('-hold_jid')
                args.append(dep)                                
            if int(array)>0:
                args.append('-t')
                args.append('1-' + str(array))                            
            args.append('-l')
            args.append('h_rt=' + time)            
                                    
        args.append(script)
        if env == 'hpc' or env == 'empty_hpc':            
            args.append(runparams['pdb'])         #1
            args.append(runparams['name'])        #2
            args.append(runparams['split'])                #3
            args.append(runparams['mutation'])    #4   
            args.append(runparams['variant'])     #5
            args.append(runparams['variantfile']) #6
            args.append(runparams['repairs'])     #7
        else:
            args.append(script)
            args.append('pdb='+runparams['pdb'])         #1
            args.append('name='+runparams['name'])        #2
            args.append('split='+runparams['split'])            #3
            args.append('mutation='+runparams['mutation'])    #4   
            args.append('variant='+runparams['variant'])     #5
            args.append('variantfile='+runparams['variantfile']) #6
            args.append('repairs='+runparams['repairs']) #7

        print(args)
        if env == 'hpc':
            #print('Running on hpc')
            process = subprocess.Popen(args=args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            result = process.communicate()
            print(result) #e.g. Your job 588483 ("foldx-aggregate") has been submitted
            results = result[0].split(' ')    
            jobid = results[2]    
            if '.' in jobid:
                results = jobid.split('.')
                jobid = results[0]
            dependencies[int(job)] = jobid
        elif env == 'python':
            #print('Running in python')
            dependencies[int(job)] = 0
            process = subprocess.Popen(args=args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)    
            result = process.communicate()
            print(result) #e.g. Your job 588483 ("foldx-aggregate") has been submitted
        else:
            #print('Not running')
            dependencies[int(job)] = 0

####################################################################################################
if __name__ == '__main__':
    import sys
    globals()['run_pipeline00'](sys.argv)