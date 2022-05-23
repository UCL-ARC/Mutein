"""
RSA: 19/5/22
This CI test script runs from the dataset level
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
def test_genesgene(inputs):
    inputs = addpath(inputs)
    import ga_1_datasetgenes as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_datasetprep(inputs):
    inputs = addpath(inputs)
    import ga_1_datasetpdbs as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)


######################################################################
### INPUTS
dataset="shearwater"
gene="ALL"
split=1000
vsplit=200

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
runs = "a"
runPPL("runs="+runs+"@dataset="+dataset+"@gene="+gene)


"""
prep = 1

if prep:
    import time
    tic = time.perf_counter()
    test_datasetprep("variant=*@dataset="+dataset+"@split="+str(split)+"@vsplit="+str(vsplit))
    toc = time.perf_counter()
    print("Prepared the dataset for",dataset,f"in {toc - tic:0.4f} seconds")

"""