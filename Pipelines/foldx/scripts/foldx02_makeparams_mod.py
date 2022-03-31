'''
-----------------------------
RSA 15/03/2022
-----------------------------
This file takes a pdb code (file must be located in the same directory in the format 1xyz.pdb)
It formats the pdb file into a paramater file suitable for foldx PositionScan
The output is in the same directory with the name
scanparams_1xyz.txt
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
def run_pipeline02(args):
    print('### FoldX make params job ###')
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

    pdb = params['pdb']
    jobname = params['name']
    chainid = params['chain']
    rows = int(params['split'])

    ##########################################################

    ##### Open the pdb file ################################    
    with open(thruput_path + pdb + '_rep.pdb') as f:
        pdbcontent = f.readlines()

    ##### Amino acid dictionary to convert between 3 and 1 codes
    aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}                 
                    
    ##### Go through each line and get the WT residue, then construct the foldx string to mutate it to everything, add to a list
    params_lst = []
    for line in pdbcontent: 
        line = line.strip()
        if line.startswith("ATOM"):
            linecontents = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
            # this is e.g.   ATOM        4481         CB            THR         A         716           
            #                          atom no        atom name     amino acid   chain   residue number
            atomtype = linecontents[2].strip()
            if atomtype == "CA":
                if linecontents[4].strip() == chainid:                
                    aaa = linecontents[3].strip()
                    if aaa in aa_dict:
                        aa = aa_dict[aaa]
                        mut = aa + linecontents[4].strip() + linecontents[5].strip() + 'a'#aa chain rid mutation = mutation string                  
                        params_lst.append(mut)
                    else:
                        print('!Error maybe?',aaa) #TODO think about this
                                    
    ##### Create a dataframe for the paramterfile in the number of chunks specified
    rows = int(rows)
    param_dic = {}
    param_dic['pdb'] = []
    param_dic['chain'] = []
    param_dic['mutation'] = []
    param_dic['row'] = []
    mut_str = ''

    total_muts = len(params_lst)
    chunk = int(total_muts/rows)
    remainder = int(total_muts%rows)
    # so until we get to the remainer we need chunk +1 on each row
    row_size = 0
    row = 0
    for i in range(len(params_lst)):
        mut = params_lst[i]
        if row_size == 0:
            param_dic['pdb'].append(pdb)
            param_dic['chain'].append(chainid)
            param_dic['mutation'].append(mut)
            row += 1
            param_dic['row'].append('row' + str(row))                    
        else:
            param_dic['mutation'][row-1] = param_dic['mutation'][row-1] + ',' + mut
        row_size += 1
        
        if row_size == chunk and row > remainder:
            row_size = 0        
        elif row_size == chunk+1:
            row_size = 0                
    print(total_muts,rows,chunk,row)
    ##### Turn the dictionary into a dataframe
    data_params = pd.DataFrame.from_dict(param_dic)
    filename = interim_path + 'params.txt'
    data_params.to_csv(filename,index=False,sep=' ',header=False)

