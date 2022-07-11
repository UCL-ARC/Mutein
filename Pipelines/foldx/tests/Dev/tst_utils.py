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

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
python_path = "/".join(dirs) + "/python"
lib_path = "/".join(dirs) + "/libs"
sys.path.append(python_path)
sys.path.append(lib_path)
import Paths
import remote
        
def test_clean():    
    #python ${script} $mode $pattern $WorkDir $DataDir $InstallDir $PipelineDir
    args = ["",
            'GENES', 
            'notch:notch1:x', 
            '/home/ucbtlcr/Scratch/workspace/', 
            '/home/rachel/UCL/github/MuteinData/', 
            '/home/rachel/UCL/github/Mutein/', 
            '/home/rachel/UCL/github/Mutein/Pipelines/foldx/']    
    remote.run_pipeline(args)


######################################
## Tests for the geneprot
test_clean()

