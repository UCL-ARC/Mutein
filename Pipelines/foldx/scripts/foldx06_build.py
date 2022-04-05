'''
-----------------------------
RSA 15/03/22
-----------------------------
Adapted from:
https://github.com/shorthouse-mrc/COVID_structure/blob/main/Alpha_variant/Foldx_variantcombinations_calculation.py
-----------------------------

This performs a mutation on a given list of amino acid positions on the structure
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
'''
import os
import pandas as pd
from shutil import copyfile
import helper as hlp

##### INPUTS #############################################
# The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
def run_pipeline06(args):
    print('### FoldX build job ###')
    print(args)
    iparams = hlp.inputparams(args)    
    pdb = ''
    if 'pdb' in iparams:
        pdb = iparams['pdb']
    cparams = hlp.configparams(pdb)
    print('CONFIG PARAMS',cparams,'\n')
    params = hlp.mergeparams(cparams,iparams)
    print(params)
    user = params['user']
    user, (foldxe, pythonexe, environment) = hlp.getenvironment(user)
    print(user, foldxe, pythonexe,environment,'\n')
    #####################################################     
    pdb = params['pdb']
    jobname = params['name']
    row = params['row']
    mutation_string = params['mutation']
    input_path, thruput_path, interim_path, output_path = hlp.get_make_paths(pdb,jobname)
    ############################################    
    # set up the files and directories
    pdbfile =pdb + '_rep' + str(params['repairs']) + '.pdb'
    mutations = []

    if mutation_string == ".":
        filename = thruput_path + 'variant_params.txt'
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
        # we have specified a mutation and a non-file row num
        mutations.append([mutation_string,'row'+ str(row)])    
    elif row[:3]=="row":
        # we have specified a mutation and row from the file
        mutations.append([mutation_string,row])    
    else:
        # we have specified a numerical row in the file
        filename = interim_path + 'params.txt'
        print('open',filename)
        with open(filename) as fr:
            paramscontent = fr.readlines()
            row = paramscontent[int(row)-1].strip()
            print(row)
            rowvals = row.split(' ')
            mutation = rowvals[2]
            row = rowvals[3]
            mutations.append([mutation,row])  

    for mut,row in mutations:
        # put mutation into a file
        row_path = interim_path + row + '/'
        print('### ... change directory',row_path)
        params['thisrow'] = row
        params['thismut'] = mut
        hlp.goto_job_dir(row_path,args,params,'_inputs06')        
        print('### ... copying file',thruput_path + pdbfile,row_path + pdbfile)
        copyfile(thruput_path + pdbfile,row_path + pdbfile)
        mut_fl = 'individual_list.txt'
        mut_log = 'buildmodel.log'
        with open(mut_fl,'w') as fw:
            fw.write(mut)        
        #~/UCL/libs/foldx5/foldx --command=BuildModel --ionStrength=0.05 --pH=7
        #  --water=CRYSTAL --vdwDesign=2 --pdbHydrogens=false --numberOfRuns=15 --mutant-file=mutations.txt
        #  --pdb=6vxx_rep.pdb > buildmodel.log
        foldxcommand = foldxe + ' --command=BuildModel'
        foldxcommand += ' --ionStrength=0.05'
        foldxcommand += ' --pH=7'
        foldxcommand += ' --water=CRYSTAL'
        foldxcommand += ' --vdwDesign=2'
        foldxcommand += ' --pdbHydrogens=false'
        foldxcommand += ' --numberOfRuns=15'
        foldxcommand += ' --mutant-file=' + mut_fl
        foldxcommand += ' --pdb=' + pdbfile
        foldxcommand += ' > ' + mut_log
        print()
        print(foldxcommand)        
        print()
        if 'empty' not in environment:    
            os.system(foldxcommand)

if __name__ == '__main__':
    import sys
    globals()['run_pipeline06'](sys.argv)

