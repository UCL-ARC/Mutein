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

ab,c,de,f,g,h = 1,0,0,0,0,0

timea = datetime.now()
if ab:
    timec = datetime.now()            
    print("### MUTEIN INSTALL TEST - PREPARING GENE ###")
    print("--- ETA= 0:00:08")
    ppl.run_pipeline(["",f"runs=a@dataset=one@install_dir={inst_d}@data_dir={data_d}"])
    ppl.run_pipeline(["",f"runs=b@dataset=one@install_dir={inst_d}@data_dir={data_d}"])
    timed = datetime.now()            
    print("Time taken = ", timed-timec)
    print("Problems with this are usually html request problems ewith bioservices")

################################################
timeb = datetime.now()            
print("Total time taken = ", timeb-timea)
            