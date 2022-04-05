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
import Arguments

def run_pipeline01(args):

    ##### INPUTS #############################################
    # The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
    print('### FoldX repair job ###')
    argus = Arguments.Arguments(args)        
    repair_path = argus.arg('interim_path') + 'repair' + argus.arg('repairs') + '/'
    argus.params['repair_path'] = repair_path
    hlp.goto_job_dir(argus.arg('repair_path'),args,argus.params,'_inputs01')    
    ############################################
    pdbfile = argus.arg('pdb') +'.pdb'            
    # Set up files (retain copy of original)
    numRepairs = int(argus.arg('repairs'))
    repairinnames = []
    repairoutnames = []
    for r in range(numRepairs+1):
        repairinnames.append(argus.arg('pdb') + '_' + str(r) + '.pdb')    
        repairoutnames.append(argus.arg('pdb') + '_' + str(r) + '_Repair.pdb')    
    repairinnames[numRepairs] = argus.arg('pdb') + '_rep' + str(numRepairs) + '.pdb'
    #### there are 2 files we need in the interim directory, pdb file rotabase, but rotabase is only needed for foldx4 and NOT needed for foldx5
    print('### ... copying file',pdbfile,argus.arg('input_path') + repairinnames[0],'... ###')
    copyfile(argus.arg('input_path') + '/' + pdbfile, argus.arg('repair_path') + repairinnames[0])        
    
    repairBaseA = argus.arg('foldxe') + ' --command=RepairPDB --pdb='
    repairBaseB = " --ionStrength=0.05 --pH=7 --vdwDesign=2 --pdbHydrogens=false"

    # Run repair number one
    for r in range(numRepairs):
        repaircommand = repairBaseA + repairinnames[r] + repairBaseB + ' > repair_' + str(r) + '.txt'
        #repaircommand = repairBaseA + repairinnames[r]
        print('### ... repair command #',r,repaircommand)
        os.system(repaircommand)    
        print('### ... copying file',argus.arg('repair_path') + repairoutnames[r],argus.arg('repair_path') + repairinnames[r+1])
        copyfile(repairoutnames[r],repairinnames[r+1])

    # copy the final repaired file to our main interim directory
    print('### ... copying file',argus.arg('repair_path') + repairoutnames[r], argus.arg('thruput_path')+repairoutnames[r])
    copyfile(argus.arg('repair_path') + repairinnames[numRepairs], argus.arg('thruput_path')+repairinnames[numRepairs])

    print('### COMPLETED FoldX repair job ###')

##########################################################################################   
if __name__ == '__main__':
    import sys
    globals()['run_pipeline01'](sys.argv)
