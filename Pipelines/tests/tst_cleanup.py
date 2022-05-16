# test_dummy.py
"""
RSA 21/4/2022
---------------------------
Any function thatstarts test_ in a file that startes test_ will be run automatically following a pull request.
This tests the outer pipeline
---------------------------

"""
import os
import sys

def addpath():        
    dirs = (os.path.dirname(os.path.realpath(__file__))).split("/")[:-2]
    newpath = "/".join(dirs)        
    print("Adding sys path=", newpath)
    sys.path.append(newpath)
        
def test_clean():
    addpath()
    # There are 7 arguments
    #install_dir = args[1]   # 1) the executable installation directory, the root directory of the peipeline    
    #working_dir = args[2]   # 1) the working dir, the root that the data output and input lives in
    
    import Pipelines.libs.pipeline_delete as pipe
    args = [
        "",
        "/home/rachel/UCL/github/Mutein/",
        "/home/rachel/UCL/github/MuteinData/",                
    ]
    pipe.run_pipeline(args)


######################################
## Tests for the geneprot
test_clean()

