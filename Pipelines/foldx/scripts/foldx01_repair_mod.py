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
from shutil import copyfile
import helper as hlp

def run_pipeline01(args):

    ##### INPUTS #############################################
    # The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
    print('### FoldX repair job ###')
    print(args)    
    # combine config and job input params
    iparams = hlp.inputparams(args)    
    pdb = ''
    if 'pdb' in iparams:
        pdb = iparams['pdb']
    cparams = hlp.configparams(pdb)    
    params = hlp.mergeparams(cparams,iparams)
    print(params)
    user = params['user']
    user, (foldxe, pythonexe, environment) = hlp.getenvironment(user)
    print(user, foldxe, pythonexe,environment)
    pdb = params['pdb']
    jobname = params['name']
    input_path, thruput_path, interim_path, output_path = hlp.get_make_paths(pdb,jobname)
    ############################################
    pdbfile = pdb +'.pdb'            
    # Set up files (retain copy of original)
    numRepairs = 5
    repairinnames = []
    repairoutnames = []
    for r in range(numRepairs+1):
        repairinnames.append(pdb + '_' + str(r) + '.pdb')    
        repairoutnames.append(pdb + '_' + str(r) + '_Repair.pdb')    
    repairinnames[numRepairs] = pdb + '_rep.pdb'

    #### there are 2 files we need in the interim directory, pdb file rotabase, but rotabase is only needed for foldx4 and NOT needed for foldx5
    print('### ... copying file',pdbfile,input_path + repairinnames[0],'... ###')
    copyfile(input_path + '/' + pdbfile, interim_path + repairinnames[0])
    # Now change into the interim directory for the work
    print('### ... changing directory to',interim_path)
    os.chdir(interim_path)

    repairBaseA = foldxe + ' --command=RepairPDB --pdb='
    repairBaseB = " --ionStrength=0.05 --pH=7 --vdwDesign=2 --pdbHydrogens=false"

    # Run repair number one
    for r in range(numRepairs):
        repaircommand = repairBaseA + repairinnames[r] + repairBaseB + ' > repair_' + str(r) + '.txt'
        #repaircommand = repairBaseA + repairinnames[r]
        print('### ... repair command #',r,repaircommand)
        os.system(repaircommand)    
        print('### ... copying file',interim_path + repairoutnames[r],interim_path + repairinnames[r+1])
        copyfile(repairoutnames[r],repairinnames[r+1])

    # copy the final repaired file to our main interim directory
    print('### ... copying file',interim_path + repairoutnames[r], thruput_path+repairoutnames[r])
    copyfile(interim_path + repairinnames[numRepairs], thruput_path+repairinnames[numRepairs])

    print('### COMPLETED FoldX repair job ###')
