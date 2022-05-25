"""
-----------------------------
RSA 24/05/22
-----------------------------

This simply views a file
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
import os
import pathlib
import pwd
import pandas as pd
from shutil import copyfile
from os.path import exists
import subprocess
from datetime import datetime

# import from the shared library in Mutein/Pipelines/shared/lib
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)

import Arguments
import Paths
from os import listdir
from os.path import isfile, join
import os


def checkResults(ddg,bm,ps):
    if exists(ddg):
        timeA = pathlib.Path(ddg).stat().st_mtime            
        print("DDG Background file was created at",datetime.fromtimestamp(timeA).strftime('%d-%m-%y-%H:%M'))
    else:
        print("!!!DDG Background file does not exist",ddg)
    if exists(bm):
        timeB = pathlib.Path(bm).stat().st_mtime
        print("Variant buildmodel file was created at",datetime.fromtimestamp(timeB).strftime('%d-%m-%y-%H:%M'))
    else:
        print("!!!Variant buildmodel file does not exist",bm)
    if exists(ps):
        timeC = pathlib.Path(ps).stat().st_mtime
        print("Variant posscan file was created at",datetime.fromtimestamp(timeC).strftime('%d-%m-%y-%H:%M'))
    else:
        print("!!!Variant posscan file does not exist",ps)



def run_pipeline(args):    
    now = datetime.now()
    current_time = now.strftime("%d-%m-%y@%H.%H.%S")
    print("Mutein remote script:",current_time)
    
    ##############################################        
    mode = args[1]
    pattern = args[2]
    WorkDir = args[3]
    DataDir = args[4]
    InstallDir = args[5]
    PipelineDir = args[6]

    print("Mode=",mode)
    print("Pattern=",pattern)
    print("")
    if mode == "GENES":
        dataset_gene_pdb=pattern.split(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset)
        filename = path.inputs + "genes_pdb_list.csv"
        if exists(filename):
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for line in lines:
                    print(line.strip())
        else:
            print("The dataset has not been prepared - no genes list",filename)
    elif mode == "PDBS":
        dataset_gene_pdb=pattern.split(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene)
        
        print("Check results files for gene")
        filenameA = path.outputs + "ddg_background.csv"
        filenameB = path.outputs + "ddg_variant_bm.csv"
        filenameC = path.outputs + "ddg_variant_ps.csv"
        checkResults(filenameA,filenameB,filenameC)

        print("\nCheck pdb list\n")        
        filename = path.outputs + "pdb_tasklist.csv"
        if exists(filename):
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for line in lines:
                    print(line.strip())
        else:
            print("The pdbs have not been prepared - no pdb list",filename)
    elif mode == "PDB":
        dataset_gene_pdb=pattern.split(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene,pdb=pdb)
        print("Check results files for pdb")
        filenameA = path.outputs + "ddg_background.csv"
        filenameB = path.outputs + "ddg_variant_bm.csv"
        filenameC = path.outputs + "ddg_variant_ps.csv"
        checkResults(filenameA,filenameB,filenameC)
        





    
    

    

##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
