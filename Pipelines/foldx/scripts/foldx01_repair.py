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
params = {}
params = hlp.inputparams(sys.argv,params)
print(params)
'''

pdb = sys.argv[1]
jobname = sys.argv[2]
realOrTest = sys.argv[3]
environment = sys.argv[4]
############################################
foldx_exe = 'foldx'
if environment == 'RSA':
    foldx_exe = '~/UCL/libs/foldx/foldx' # on my laptop RSA

pdbfile = pdb +'.pdb'

# we want to work in the node directory first, the main pdb input file is a 1-off and lives in github (at the moment)
dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
print('### ... changing directory to',dir_path)
os.chdir(dir_path)
results_dir = dir_path + 'results/' + jobname + '/' 
resultrepair_dir = results_dir + 'repair/'

# Set up the directory and any files we need (may not be in the path)
if not os.path.exists(results_dir):
    os.mkdir(results_dir)
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

#### there are 2 files we need in the results directory, pdb file rotabase
print('### ... copying file',pdbfile,results_dir + repairinnames[0],'... ###')
copyfile(pdbfile, resultrepair_dir + repairinnames[0])
#print('### ... copying file','rotabase.txt',results_dir + 'rotabase.txt')
#copyfile('rotabase.txt', resultrepair_dir + 'rotabase.txt')
# Now change into the results directory for the work
print('### ... changing directory to',results_dir)
os.chdir(resultrepair_dir)

repairBaseA = foldx_exe + ' --command=RepairPDB --pdb='
repairBaseB = " --ionStrength=0.05 --pH=7 --vdwDesign=2 --pdbHydrogens=false"

# Run repair number one
for r in range(numRepairs):
    repaircommand = repairBaseA + repairinnames[r] + repairBaseB + ' > repair_' + str(r) + '.txt'
    #repaircommand = repairBaseA + repairinnames[r]
    print('### ... repair command #',r,repaircommand)
    os.system(repaircommand)    
    print('### ... copying file',results_dir + repairoutnames[r],results_dir + repairinnames[r+1])
    copyfile(repairoutnames[r],repairinnames[r+1])

# copy the final repaired file to our main results directory
print('### ... copying file',resultrepair_dir + repairoutnames[r], dir_path+repairoutnames[r])
copyfile(resultrepair_dir + repairinnames[numRepairs], results_dir+repairinnames[numRepairs])

print('### COMPLETED FoldX repair job ###')
'''