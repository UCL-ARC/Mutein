'''
------------------------
RSA 29/03/22
------------------------
Pipeline script for foldx job on Myriad
------------------------
'''

import sys
import os
import subprocess
import helper as hlp

print('#### FOLDX PIPELINE - batch creation ####')
# Make sure out current directory is where the script lives
dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
print('## ... changing directory to',dir_path)
os.chdir(dir_path)
##### INPUTS #############################################
# 1= repairing pdb
# 2= making param file for splits
# 3= performing the position scan ddg mutations (parallel)
# 4= aggregating 3
# 5= making params file for variants
# 6 = performing variant ddg (parallel)
# 7 = aggregating 6
jobparams = hlp.inputparams(sys.argv)
cfgparams = hlp.configparams(jobparams['configfile'])
# the job overrides the config
runparams = hlp.mergeparams(cfgparams, jobparams)
user = runparams['user']
user, (foldxe, pythonexe, env) = hlp.getenvironment(user)
print(user, foldxe, pythonexe,env)

###### For the purposes of testing override params as required #####
runparams['jobs'] = '1'
print('final params',runparams)

# script extension is either sh or py depending on bash or python environment
ext = '.sh'
if env == 'python':
    ext = '.py'

#############################################################
dependencies = {1:-1,2:-1,3:-1,4:-1,5:-1,6:-1,7:-1}
runs = []
if '1' in runparams['jobs']:
    runs.append([1,'qsub',"./foldx01_repair" + ext,-1])
if '2' in runparams['jobs']:
    if '1' in runparams['jobs']:
        runs.append([2,'qsub',"./Sh02_Myriad_makeparams" + ext,1])
    else:
        runs.append([2,'qsub',"./Sh02_Myriad_makeparams" + ext,-1])

override_array = False
if '3' in runparams['jobs']:
    dep = -1
    if '2' in runparams['jobs']:
        dep = 2
    if runparams['mutation'] == '.':
        override_array = True
        runs.append([3,'qsub',"./Sh03_Myriad_posscan" + ext,dep])        
    else:
        runs.append([3,'qsub',"./Sh03_Myriad_posscan_one" + ext,dep])

if '4' in runparams['jobs']:    
    if '3' in runparams['jobs']:
        runs.append([4,'qsub',"./Sh04_Myriad_aggregate" + ext,3])
    else:
        runs.append([4,'qsub',"./Sh04_Myriad_aggregate" + ext,-1])

if '5' in runparams['jobs']:
    if '1' in runparams['jobs']:
        runs.append([5,'qsub',"./Sh05_Myriad_makevparams" + ext,1])
    else:
        runs.append([5,'qsub',"./Sh05_Myriad_makevparams" + ext,-1])

if '6' in runparams['jobs']:
    if '5' in runparams['jobs']:
        runs.append([6,'qsub',"./Sh06_Myriad_build" + ext,5])
    else:
        runs.append([6,'qsub',"./Sh06_Myriad_build" + ext,-1])

if '7' in runparams['jobs']:
    if '6' in runparams['jobs']:
        runs.append([7,'qsub',"./Sh07_Myriad_vaggregate" + ext,6])
    else:
        runs.append([7,'qsub',"./Sh07_Myriad_vaggregate" + ext,-1])

for job,exe,script,dependency in runs:
    print(job,exe,script,dependency)
    if env == 'hpc':
        os.system('chmod +x ' + script)
    args = []
    if env == 'python':
        args.append(pythonexe)        
    else:
        args.append('qsub')
        
        if dependency > -1:
            dep = dependencies[int(dependency)]
            args.append('-hold_jid')
            args.append(dep)        
        
        if override_array and job == 3:        
            args.append('-t')
            args.append('1-' + str(runparams['split']))    
        
        if job == 3 and runparams['time'] != '.':
            args.append('-l')
            args.append('h_rt=' + runparams['time'])    
        
        if job == 6:        
            args.append('-t')
            args.append('1-' + str(runparams['combos']))
                    
    args.append(script)
    args.append(runparams['pdb'])         #1
    args.append(runparams['name'])        #2
    args.append(runparams['split'])       #3
    args.append(runparams['mutation'])    #4   
    args.append(runparams['variant'])     #5
    args.append(runparams['variantfile']) #6


    print(args)
    if env == 'hpc':
        print('Running on hpc')
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
        print('Running in python')
        dependencies[int(job)] = 0
        process = subprocess.Popen(args=args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)    
        result = process.communicate()
        print(result) #e.g. Your job 588483 ("foldx-aggregate") has been submitted
    else:
        print('Not running')
        dependencies[int(job)] = 0
    