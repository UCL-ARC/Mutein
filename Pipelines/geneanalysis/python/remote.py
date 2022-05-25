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

    if mode == "GENES":
        dataset_gene_pdb=pattern.spit(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene,pdb=pdb)
        filename = path.inputs + "genes_pdb_list.csv"
        if exists(filename):
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for line in lines:
                    print(line)
        else:
            print("The dataset has not been prepared - no genes list")



    
    

    

##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
