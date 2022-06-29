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

b,c,de,f,g,h = 0,0,0,1,0,0

timea = datetime.now()

if b:
    timec = datetime.now()            
    print("### MUTEIN INSTALL TEST - DOWNLOAD FROM WEB ###")
    print("--- ETA= 0:00:03")
    ppl.run_pipeline(["",f"runs=b@pdb=1pb5@install_dir={inst_d}@data_dir={data_d}"])
    timed = datetime.now()            
    print("Time taken = ", timed-timec)
if c:
    timec = datetime.now()            
    print("### MUTEIN INSTALL TEST - REPAIRING PROTEIN STRUCTURE ###")
    print("--- ETA= 0:00:08")
    ppl.run_pipeline(["",f"runs=c@pdb=1pb5@repairs=1@repair_from=0@install_dir={inst_d}@data_dir={data_d}"])
    timed = datetime.now()            
    print("Time taken = ", timed-timec)
if de:
    timec = datetime.now()            
    print("### MUTEIN INSTALL TEST - SPLITTING FOR PARALLELISATION ###")
    print("--- ETA= 0:00:00")
    ppl.run_pipeline(["",f"runs=d@pdb=1pb5@repairs=x@chunk=500@variant=*@vchunk=10@install_dir={inst_d}@data_dir={data_d}"])
    ppl.run_pipeline(["",f"runs=e@pdb=1pb5@repairs=x@chunk=10@variant=*@vchunk=10@install_dir={inst_d}@data_dir={data_d}"])
    timed = datetime.now()            
    print("Time taken = ", timed-timec)
if f:
    timec = datetime.now()            
    print("### MUTEIN INSTALL TEST - RUNNING BACKGROUND TASKS ###")
    print("--- ETA= 0:16:00")
    ppl.run_pipeline(["",f"runs=f@pdb=1pb5@repairs=x@task=1@install_dir={inst_d}@data_dir={data_d}"])
    ppl.run_pipeline(["",f"runs=f@pdb=1pb5@repairs=x@task=2@install_dir={inst_d}@data_dir={data_d}"])
    timed = datetime.now()            
    print("Time taken = ", timed-timec)
if g:
    timec = datetime.now()            
    print("### MUTEIN INSTALL TEST - RUNNING VARIANT TASKS ###")
    print("--- ETA= 0:00:20")
    ppl.run_pipeline(["",f"runs=g@pdb=1pb5@repairs=x@task=1@variant=*@vchunk=10@install_dir={inst_d}@data_dir={data_d}"])
    ppl.run_pipeline(["",f"runs=g@pdb=1pb5@repairs=x@task=2@variant=*@vchunk=10@install_dir={inst_d}@data_dir={data_d}"])
    timed = datetime.now()            
    print("Time taken = ", timed-timec)
if h:
    timec = datetime.now()            
    print("### MUTEIN INSTALL TEST - AGGREGATING ###")
    print("--- ETA= 0:00:01")
    ppl.run_pipeline(["",f"runs=h@pdb=1pb5@repairs=x@install_dir={inst_d}@data_dir={data_d}"])
    timed = datetime.now()            
    print("Time taken = ", timed-timec)



################################################
timeb = datetime.now()            
print("Total time taken = ", timeb-timea)
            