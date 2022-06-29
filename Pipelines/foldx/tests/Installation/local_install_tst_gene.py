"""
RSA: 19/5/22
This CI test script runs from the pdb level
It enables debugging of the scripts as if run from a batch

"""
import sys, os
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

######## INPUTS #####################################
envs = env.getenvironment() # WARNING YOU NEED TO HAVE ADDED YOURSELF HERE
inst_d = envs["install_dir"]
data_d = envs["data_dir"]
dataset="mouse2"

a,b = 0,1

### Note the HTML links do not work if you have CISCO connect turned on.###

timea = datetime.now()
if a:
    timec = datetime.now()            
    print("### MUTEIN INSTALL TEST - PREPARING GENE ###")
    print("--- ETA= 0:01:00")
    #ppl.run_pipeline(["",f"runs=a@dataset=cutdown@install_dir={inst_d}@data_dir={data_d}"])
    ppl.run_pipeline(["",f"runs=a@dataset={dataset}@install_dir={inst_d}@data_dir={data_d}"])
    timed = datetime.now()            
    print("Time taken = ", timed-timec)
    print("Problems with this are usually html request problems ewith bioservices")

if b:
    timec = datetime.now()            
    print("### MUTEIN INSTALL TEST - PREPARING PDBS ###")
    print("--- ETA= 0:00:08")
    #ppl.run_pipeline(["",f"runs=a@dataset=cutdown@install_dir={inst_d}@data_dir={data_d}"])
    ppl.run_pipeline(["",f"runs=b@dataset={dataset}@install_dir={inst_d}@data_dir={data_d}"])
    timed = datetime.now()            
    print("Time taken = ", timed-timec)
    
################################################
timeb = datetime.now()            
print("Total time taken = ", timeb-timea)
            