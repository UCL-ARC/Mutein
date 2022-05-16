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
import pandas as pd
from shutil import copyfile
from os.path import exists

# import from the shared library in Mutein/Pipelines/shared/lib
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)
import Paths
import Arguments
import Config
import Analysis
import FileDf


def run_pipeline(args):
    print("### Delete unnecessary log files ###")
    print(args)
    ##############################################
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")
    sys.path.append(install_dir)
    sys.path.append(install_dir + "/Pipelines")
    sys.path.append(install_dir + "/Pipelines/libs")
    data_dir = argus.arg("data_dir")
    install_dir = argus.arg("install_dir")    
    errorfile = argus.arg("errorfile")
    outfile = argus.arg("outfile")
    print("data_dir",data_dir)
    print("install_dir",install_dir)
    print("errorfile",errorfile)
    print("outfile",outfile)
    
    
    #work_path = gene_path.gene_outputs + "agg/"
    #argus.params["work_path"] = work_path
    #gene_path.goto_job_dir(argus.arg("work_path"), args, argus.params, "_inputs01")
    ############################################
    
                            
    print("### COMPLETED GENE STITCH job ###")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
