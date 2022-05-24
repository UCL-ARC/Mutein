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
    print("### Delete unnecessary log files ###")
    print(args)
    ##############################################    
    argus = Arguments.Arguments(args)    
    jobid = argus.arg("JOBID")
    homeuser = pwd.getpwuid(os.getuid())[0]
    scratch_dir = "/home/" + homeuser + "/Scratch/workspace/"    
    print("scratch_dir",scratch_dir)
    
    onlyfiles = [f for f in listdir(scratch_dir) if isfile(join(scratch_dir, f))]        
    for file in onlyfiles:
        filename =scratch_dir+file 
        if jobid in filename:
            if exists(filename):                     
                with open(filename) as fr:
                    lines = fr.readlines()
                    for line in lines:                
                        print(line)                
            


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
