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
    import ga_genestogene as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_datasetprep(inputs):
    inputs = addpath(inputs)
    import ga_1_0_datasetprep as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)


######################################################################
### INPUTS
dataset="shearwater"
split=1000
vsplit=200

#test_genesgene("dataset="+dataset)

prep = 1

if prep:
    import time
    tic = time.perf_counter()
    test_datasetprep("variant=*@dataset="+dataset+"@split="+str(split)+"@vsplit="+str(vsplit))
    toc = time.perf_counter()
    print("Prepared the dataset for",dataset,f"in {toc - tic:0.4f} seconds")

