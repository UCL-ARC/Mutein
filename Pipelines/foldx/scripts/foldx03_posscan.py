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
import Arguments

##### INPUTS #############################################
# The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)

def run_pipeline03(args):    
    print('### FoldX position scan job ###')
    print(args)    
    argus = Arguments.Arguments(args)            
    
    pdb = argus.arg('pdb')        
    row = argus.arg('row')
    mutation_string = argus.arg('mutation')            
    ############################################
    pdbfile =pdb + '_rep' + str(argus.arg('repairs')) + '.pdb'            
    mutations = []
    # row=. means all, row=1:n means an explicit row, row=0 means the mutation string has been passd in explicitly
    if mutation_string == ".":
        filename = argus.arg('interim_path') + 'params.txt'
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
        filename = argus.arg('interim_path') + 'params.txt'
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
        print(mut,row)

        row_path = argus.arg('interim_path') + row + '/'
        print('### ... change directory',row_path)
        argus.params['thisrow'] = row
        argus.params['thismut'] = mut
        hlp.goto_job_dir(row_path,args,argus.params,'_inputs03')        
        print('### ... copying file',argus.arg('thruput_path') + pdbfile,row_path + pdbfile)
        copyfile(argus.arg('thruput_path') + pdbfile,row_path + pdbfile)
        
        
        foldxcommand = argus.arg('foldxe') + ' --command=PositionScan'
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
        print("RealOrTest=",argus.arg('environment'))
        if argus.arg('environment') != 'inputs':    
            os.system(foldxcommand)

##########################################################################################   
if __name__ == '__main__':
    import sys
    globals()['run_pipeline03'](sys.argv)