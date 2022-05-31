"""
RSA: 19/5/22
This CI test script runs from the pdb level
It enables debugging of the scripts as if run from a batch

"""
import sys, os
# import from the shared library in Mutein/Pipelines/shared/lib
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)
retpath = "/".join(dirs) + "/geneanalysis/python"
sys.path.append(retpath)
import Paths

def addpath(inputs):        
    print(sys.path)
    inputs += "@install_dir=/home/rachel/UCL/github/Mutein/"
    inputs += "@data_dir=/home/rachel/UCL/github/MuteinData/"
    return inputs

######################################################################
def test_proteinprep(inputs):
    inputs = addpath(inputs)
    import ga_3_0_pdbprep as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_proteinsplit(inputs):
    inputs = addpath(inputs)
    import ga_3_pdbbackparams as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_proteinvsplit(inputs):
    inputs = addpath(inputs)
    import ga_3_pdbvarparams as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_proteinrepair(inputs):
    inputs = addpath(inputs)
    import ga_3_proteinrepair as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_proteintasks(inputs):
    inputs = addpath(inputs)
    import ga04a_posscan as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_proteinvtasks(inputs):
    inputs = addpath(inputs)
    import ga04b_singlescan as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_proteinagg(inputs):
    inputs = addpath(inputs)
    import ga05a_aggddg as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_proteinvagg(inputs):
    inputs = addpath(inputs)
    import ga05b_singlesagg as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_proteindblagg(inputs):
    inputs = addpath(inputs)
    import ga_2_genestitch as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)
######################################################################
### INPUTS
dataset=""
gene=""
pdb="1pb5"

repairs=2
split=2
vsplit=2
task=1
runs = "hi"
""" 
    if "a" in runs:        print("Mutein: Preparing genes")        
    if "b" in runs:        print("Mutein: Preparing pdbs")        
    if "c" in runs:        print("Mutein: Repairing pdbs")        
    if "d" in runs:        print("Mutein: Making background param file")        
    if "e" in runs:        print("Mutein: Making variant param file")        
    if "f" in runs:        print("Mutein: Running background tasks")        
    if "g" in runs:        print("Mutein: Running variant tasks")        
    if "h" in runs:        print("Mutein: Aggregating background tasks")        
    if "i" in runs:        print("Mutein: Aggregating variant tasks")        
    if "j" in runs:        print("Mutein: Aggregating gene tasks")        
"""

# which steps of the pipeline to run
def runPPL(inputs):
    inputs = addpath(inputs)
    import ga__runner as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)
runPPL("runs="+runs+"@dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@repairs="+str(repairs)+"@task=" +str(task)+"@chunk="+str(split)+"@variant=*@vchunk="+str(vsplit))

