'''
-----------------------------
RSA 15/03/22
-----------------------------
Adapted from:
https://github.com/shorthouse-mrc/COVID_structure/blob/main/Foldx_positionscanall.py
-----------------------------

This performs a mutation on a given list of amino acid positions on the structure
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
'''
#import sys
import os
import pandas as pd
from shutil import copyfile
import helper as hlp

##### INPUTS #############################################
# The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)

def run_pipeline03(args):    
    print('### FoldX position scan job ###')
    print(args)    
    # combine config and job input params
    # additional inputs special to this job are row and mutation
    ##############################################
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
    row = params['row']
    mutation_string = params['mutation']
    input_path, thruput_path, interim_path, output_path = hlp.get_make_paths(pdb,jobname)
    ############################################
    pdbfile = pdb +'_rep.pdb'
    print('### ... change directory',thruput_path)                     
    mutations = []
    # row=. means all, row=1:n means an explicit row, row=0 means the mutation string has been passd in explicitly
    if row == ".":
        filename = interim_path + 'params.txt'
        print('open',filename)
        with open(filename) as fr:
            paramscontent = fr.readlines()        
            for row in paramscontent:
                row = row.strip()
                print(row)
                rowvals = row.split(' ')
                mutation = rowvals[2]
                row = rowvals[3]
                mutations.append([mutation,row])
    elif row[0]=="0":
        mutations.append([mutation_string,'row'+ str(row)])    
    else:
        filename = interim_path + 'params.txt'
        print('open',filename)
        with open(filename) as fr:
            paramscontent = fr.readlines()
            row = paramscontent[int(row-1)].strip()
            print(row)
            rowvals = row.split(' ')
            mutation = rowvals[2]
            row = rowvals[3]
            mutations.append([mutation,row])  

    

    for mut,row in mutations:
        print(mut,row)

        row_path = interim_path + row + '/'
        hlp.goto_job_dir(row_path,args,params,'_inputs03')        
        print('### ... copying file',thruput_path + pdbfile,row_path + pdbfile)
        copyfile(thruput_path + pdbfile,row_path + pdbfile)
        print('### ... change directory',row_path)
        os.chdir(row_path)

        foldxcommand = foldxe + ' --command=PositionScan'
        foldxcommand += ' --ionStrength=0.05'
        foldxcommand += ' --pH=7'
        foldxcommand += ' --water=CRYSTAL'
        foldxcommand += ' --vdwDesign=2'
        foldxcommand += ' --pdbHydrogens=false'
        #foldxcommand += ' --output-dir=' + results_dir
        #foldxcommand += ' --pdb-dir=' + results_dir
        foldxcommand += ' --pdb=' + pdbfile
        foldxcommand += ' --positions=' + mut


        print(foldxcommand)
        print("RealOrTest=",environment)
        if environment != 'empty':    
            os.system(foldxcommand)

##########################################################################################   
if __name__ == '__main__':
    import sys
    globals()['run_pipeline03'](sys.argv)