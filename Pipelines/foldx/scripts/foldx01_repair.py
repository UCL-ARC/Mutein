'''
-----------------------------
RSA 15/03/2022
-----------------------------
Adapted from
https://github.com/shorthouse-mrc/COVID_structure/blob/main/Foldx_repair6.py
-----------------------------

This file takes a pdb code (file must be located in the same directory in the format 1xyz.pdb)
It formats the pdb file into a paramater file suitable for foldx PositionScan
The output is in the same directory with the name
scanparams_1xyz.txt
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
'''
import os
import sys
from shutil import copyfile
import helper as hlp

##### INPUTS #############################################
# The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
print('### FoldX repair job ###')
print(sys.argv)
iparams = hlp.inputparams(sys.argv)
cparams = hlp.configparams('')
params = hlp.mergeparams(cparams,iparams)
print(params)
user = params['user']
user, (foldxe, pythonexe, environment) = hlp.getenvironment(user)
print(user, foldxe, pythonexe,environment)

pdb = params['pdb']
jobname = params['name']
############################################

pdbfile = pdb +'.pdb'

# we want to work in the node directory first, the main pdb input file is a 1-off and lives in github (at the moment)

#\\wsl$\Ubuntu\home\rachel\UCL\github\Mutein\Pipelines\foldx\scripts
#\\wsl$\Ubuntu\home\rachel\UCL\github\Mutein\Pipelines\foldx\inputs\6vxx

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path = dir_path[:-7]
pdb_path = dir_path + 'inputs/'+pdb + '/'
results_path = dir_path + 'results/'+jobname + '/'
resultrepair_dir = results_path + 'repair/'
print('script path',dir_path)
print('pdb path',pdb_path)
print('results path',results_path)
print('repair path',resultrepair_dir)
# All the files MUST exist already apart from the jobname folders - results_path and resultrepair_path
# Set up the directory and any files we need (may not be in the path)
if not os.path.exists(dir_path + 'results/'):
    os.mkdir(dir_path + 'results/')
if not os.path.exists(dir_path + 'thruputs/'):
    os.mkdir(dir_path + 'thruputs/')
if not os.path.exists(dir_path + 'thruputs/' + pdb + '/'):
    os.mkdir(dir_path + 'thruputs/' + pdb + '/')
os.chdir(dir_path + 'results/')
if not os.path.exists(results_path):
    os.mkdir(results_path)
if not os.path.exists(resultrepair_dir):
    os.mkdir(resultrepair_dir)
# Set up files (retain copy of original)
numRepairs = 5
repairinnames = []
repairoutnames = []
for r in range(numRepairs+1):
    repairinnames.append(pdb + '_' + str(r) + '.pdb')    
    repairoutnames.append(pdb + '_' + str(r) + '_Repair.pdb')    
repairinnames[numRepairs] = pdb + '_rep.pdb'

#### there are 2 files we need in the results directory, pdb file rotabase, but rotabase is only needed for foldx4 and NOT needed for foldx5
print('### ... copying file',pdbfile,pdb_path + repairinnames[0],'... ###')
copyfile(pdb_path + '/' + pdbfile, resultrepair_dir + repairinnames[0])
# Now change into the results directory for the work
print('### ... changing directory to',resultrepair_dir)
os.chdir(resultrepair_dir)

repairBaseA = foldxe + ' --command=RepairPDB --pdb='
repairBaseB = " --ionStrength=0.05 --pH=7 --vdwDesign=2 --pdbHydrogens=false"

# Run repair number one
for r in range(numRepairs):
    repaircommand = repairBaseA + repairinnames[r] + repairBaseB + ' > repair_' + str(r) + '.txt'
    #repaircommand = repairBaseA + repairinnames[r]
    print('### ... repair command #',r,repaircommand)
    os.system(repaircommand)    
    print('### ... copying file',resultrepair_dir + repairoutnames[r],resultrepair_dir + repairinnames[r+1])
    copyfile(repairoutnames[r],repairinnames[r+1])

# copy the final repaired file to our main results directory
print('### ... copying file',resultrepair_dir + repairoutnames[r], dir_path+repairoutnames[r])
copyfile(resultrepair_dir + repairinnames[numRepairs], dir_path + 'thruputs/' + pdb + '/'+repairinnames[numRepairs])

print('### COMPLETED FoldX repair job ###')
