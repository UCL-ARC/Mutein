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
import ga__runner as ppl
import environments as env
import remote as rem

######## INPUTS #####################################
envs = env.getenvironment() # WARNING YOU NEED TO HAVE ADDED YOURSELF HERE
inst_d = envs["install_dir"]
data_d = envs["data_dir"]

cde,untasks,fg = 0,1,0

### Note the HTML links do not work if you have CISCO connect turned on.###
# we first need to make sure that the directories are empty so the tasks have not been completed

timea = datetime.now()
if cde:
    timec = datetime.now()            
    print("### MUTEIN TEST - DOWNLOAD FROM WEB ###")
    print("--- ETA= 0:00:03")       
    ppl.run_pipeline(["",f"runs=c@dataset=cutdown@gene=cr2@repairs=0@task=1@repair_from=0@install_dir={inst_d}@data_dir={data_d}"])
    ppl.run_pipeline(["",f"runs=d@dataset=cutdown@gene=cr2@repairs=x@chunk=500@variant=*@vchunk=10@install_dir={inst_d}@data_dir={data_d}"])
    ppl.run_pipeline(["",f"runs=e@dataset=cutdown@gene=cr2@repairs=x@chunk=10@variant=*@vchunk=10@install_dir={inst_d}@data_dir={data_d}"])
    timed = datetime.now()            
    print("Time taken = ", timed-timec)

if untasks:
    timec = datetime.now()            
    print("### MUTEIN TEST - CREATING UNTASKS ###")
    print("--- ETA= 0:16:00")
    rem.run_pipeline(["",'GENEINCOMPLETE', 'cutdown:cr2:', data_d + "/log/", data_d, inst_d, inst_d+"/Pipelines/foldx/"])

if fg:
    timec = datetime.now()            
    print("### MUTEIN TEST - RUNNING UNTASKS ###")
    print("--- ETA= 0:16:00")
    ppl.run_pipeline(["",f"runs=f@dataset=cutdown@gene=cr2@missing=Y@repairs=x@task=1@install_dir={inst_d}@data_dir={data_d}"])
    ppl.run_pipeline(["",f"runs=g@dataset=cutdown@gene=cr2@missing=Y@repairs=x@task=1@variant=*@vchunk=10@install_dir={inst_d}@data_dir={data_d}"])
    timed = datetime.now()            
    print("Time taken = ", timed-timec)
    
################################################
timeb = datetime.now()            
print("Total time taken = ", timeb-timea)
            