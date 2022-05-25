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

# import from the shared library in Mutein/Pipelines/shared/lib
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)
import Arguments
from os import listdir
from os.path import isfile, join
import os


def run_pipeline(args):    
    print(args)
    ##############################################        
    mode = args[1]
    pattern = args[2]
    WorkDir = args[3]
    DataDir = args[4]
    InstallDir = args[5]
    PipelineDir = args[6]

    print(mode)
    print(pattern)
    print(WorkDir)
    print(DataDir)
    print(InstallDir)
    print(PipelineDir)

    

##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
