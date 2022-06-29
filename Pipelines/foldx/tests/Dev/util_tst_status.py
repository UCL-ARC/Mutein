"""
RSA: 19/5/22
This CI test script runs from the pdb level
It enables debugging of the scripts as if run from a batch

"""
import sys, os
from os import listdir
from os.path import isfile, join
from datetime import datetime
# import from the shared library in Mutein/Pipelines/shared/lib
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
python_path = "/".join(dirs) + "/python"
lib_path = "/".join(dirs) + "/libs"
sys.path.append(python_path)
sys.path.append(lib_path)
import Paths
import environments as env
import remote as rem

######## INPUTS #####################################
envs = env.getenvironment() # WARNING YOU NEED TO HAVE ADDED YOURSELF HERE
inst_d = envs["install_dir"]
data_d = envs["data_dir"]

timea = datetime.now()


timec = datetime.now()            
print("### MUTEIN TEST - CREATING UNTASKS ###")
print("--- ETA= 0:01:00")
rem.run_pipeline(["",'GENES', 'cutdown:cr2:', data_d + "/log/", data_d, inst_d, inst_d+"/Pipelines/foldx/"])
    
################################################
timeb = datetime.now()            
print("Total time taken = ", timeb-timea)
            