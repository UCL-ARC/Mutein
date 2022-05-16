"""
-----------------------------
RSA 16/05/22
-----------------------------

This aggregates the outputs from split positionscans into 1 file
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
import os
import pwd
import pandas as pd
from shutil import copyfile
from os.path import exists

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
    homeuser = pwd.getpwuid(os.getuid())[0]
    scratch_dir = "/home/" + homeuser + "/Scratch/workspace/"    
    print("scratch_dir",scratch_dir)
    
    onlyfiles = [f for f in listdir(scratch_dir) if isfile(join(scratch_dir, f))]    
    file_numbers = {}
    for file in onlyfiles:
        name = file.split(".")[0]
        number = int(file.split(".")[1][1:])
        if number not in file_numbers:
            file_numbers[number] = name

    for number,name in file_numbers.items():
        error_file = scratch_dir+name +".e" + str(number)
        out_file = scratch_dir+name +".o" + str(number)
        if exists(out_file) and exists(error_file):        
            with open(error_file) as fr:
                lines_err = fr.readlines()
            with open(out_file) as fr:
                lines_out = fr.readlines()
            if len(lines_err) == 0:
                if len(lines_out)>0:
                    if lines_out[-1] == "MUTEIN SCRIPT ENDED":                
                        os.remove(error_file)                
                        os.remove(out_file)
                        print("...removing",number,name)
            else:
                print("Errors",number,name,lines_err)                
        else:
            print("Missing files",number,name)
                    
    
    


    
    
    
    
                            
    print("### COMPLETED CLEAN UP job ###")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
