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
    # PATHS for the batch
    dir_path = os.path.dirname(os.path.realpath(__file__))
    dir_path = dir_path[:-7]
    input_path = dir_path + 'inputs/'
    
    print('script path',dir_path)
    
    # combine config and job input params
    iparams = hlp.inputparams(args)
    configfile = ''
    if 'pdb' in iparams:
        configfile = input_path + iparams['pdb'] + '/config.cfg'
    cparams = hlp.configparams(configfile)
    print('cfg file',cparams)
    params = hlp.mergeparams(cparams,iparams)
    print(params)
    user = params['user']
    user, (foldxe, pythonexe, environment) = hlp.getenvironment(user)
    print(user, foldxe, pythonexe,environment)

    pdb = params['pdb']
    jobname = params['name']
    ############################################
    pdb_path = dir_path + 'inputs/'+pdb + '/'
    interim = dir_path + 'interim/'+jobname + '/'
    resultrepair_dir = interim + 'repair/'
    print('pdb path',pdb_path)
    print('interim path',interim)
    print('repair path',resultrepair_dir)

    pdbfile = pdb +'.pdb'

    # we want to work in the node directory first, the main pdb input file is a 1-off and lives in github (at the moment)

    #\\wsl$\Ubuntu\home\rachel\UCL\github\Mutein\Pipelines\foldx\scripts
    #\\wsl$\Ubuntu\home\rachel\UCL\github\Mutein\Pipelines\foldx\inputs\6vxx

    
    # All the files MUST exist already apart from the jobname folders - results_path and resultrepair_path
    # Set up the directory and any files we need (may not be in the path)
    if not os.path.exists(dir_path + 'interim/'):
        os.mkdir(dir_path + 'interim/')
    if not os.path.exists(dir_path + 'thruputs/'):
        os.mkdir(dir_path + 'thruputs/')
    if not os.path.exists(dir_path + 'thruputs/' + pdb + '/'):
        os.mkdir(dir_path + 'thruputs/' + pdb + '/')
    os.chdir(dir_path + 'interim/')
    if not os.path.exists(interim):
        os.mkdir(interim)
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

    #### there are 2 files we need in the interim directory, pdb file rotabase, but rotabase is only needed for foldx4 and NOT needed for foldx5
    print('### ... copying file',pdbfile,pdb_path + repairinnames[0],'... ###')
    copyfile(pdb_path + '/' + pdbfile, resultrepair_dir + repairinnames[0])
    # Now change into the interim directory for the work
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

    # copy the final repaired file to our main interim directory
    print('### ... copying file',resultrepair_dir + repairoutnames[r], dir_path+repairoutnames[r])
    copyfile(resultrepair_dir + repairinnames[numRepairs], dir_path + 'thruputs/' + pdb + '/'+repairinnames[numRepairs])

    print('### COMPLETED FoldX repair job ###')
