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
    mode = argus.arg("MODE")
    homeuser = pwd.getpwuid(os.getuid())[0]
    scratch_dir = "/home/" + homeuser + "/Scratch/workspace/"    
    print("scratch_dir",scratch_dir)
    
    onlyfiles = [f for f in listdir(scratch_dir) if isfile(join(scratch_dir, f))]    
    if mode == "ALL":
        for file in onlyfiles:
            if exists(scratch_dir+file):     
                os.remove(scratch_dir + file)   

    elif mode == "CLEAN":
        file_numbers = {}
        for file in onlyfiles:
            name = file.split(".")[0]
            number = file.split(".")[1][1:]
            if len(file.split(".")) > 2:
                number+= "." + file.split(".")[2]
            
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
                        if lines_out[-1].strip() == "MUTEIN SCRIPT ENDED":                
                            os.remove(error_file)                
                            os.remove(out_file)
                            print("...removing",number,name)
                else:
                    print("Errors",number,name)                
            else:
                print("Missing files",number,name)
    
    elif mode == "RERUN":
        file_numbers = {}
        for file in onlyfiles:
            name = file.split(".")[0]
            number = file.split(".")[1][1:]
            if len(file.split(".")) > 2:
                number+= "." + file.split(".")[2]
            
            if number not in file_numbers:
                file_numbers[number] = name

        for number,name in file_numbers.items():
            error_file = scratch_dir+name +".e" + str(number)
            out_file = scratch_dir+name +".o" + str(number)
            if exists(out_file):                        
                with open(out_file) as fr:
                    lines_out = fr.readlines()
                if len(lines_out) > 1:
                    cmd = lines_out[1].strip()
                    if len(cmd)>3:
                        print("Rerunning cmd:",cmd)
                        # but if it is an array job we only want to run the single task
                        args = cmd.split(" ")
                        if "pipeline_array" in cmd:
                            cmd = cmd.replace("pipeline_array","pipeline_single")
                            task=number.split(".")[1]
                            args[-3] += "@task=" + task
                            args = args[2:]
                            args[0] = "qsub"

                        
                        process = subprocess.Popen(args=args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                        result = process.communicate()
                        print(result)                    
                    else:
                        print("Can't re-run as manually started:",out_file)
                                                                                                                                            
    print("### COMPLETED CLEAN UP job ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
