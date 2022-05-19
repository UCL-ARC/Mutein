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
def test_geneproteins(inputs):
    inputs = addpath(inputs)
    import ga_2_genetoproteins as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_generepair(inputs):
    inputs = addpath(inputs)
    import ga_2_generepair as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_makeparams(inputs):
    inputs = addpath(inputs)
    import ga_2_genebackparams as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_makevparams(inputs):
    inputs = addpath(inputs)
    import ga_2_genevarparams as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)


######################################################################
### INPUTS
dataset=""
gene="notch1"
#test_geneproteins("dataset="+dataset+"@gene="+gene)

repairs=1
task=1
#test_generepair("dataset="+dataset+"@gene="+gene+"@repairs="+str(repairs)+"@task="+str(task))

split=100
test_makeparams("dataset="+dataset+"@gene="+gene+"@split="+str(split))

vsplit=20
test_makevparams("dataset="+dataset+"@gene="+gene+"@vsplit="+str(vsplit))
