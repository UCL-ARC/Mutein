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
def test_proteinrepair(inputs):
    inputs = addpath(inputs)
    import ga_3_proteinrepair as ppl
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


######################################################################
### INPUTS
dataset=""
gene="notch1"
pdb="1pb5"

repairs=1
#test_proteinrepair("dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@repairs="+str(repairs))

split=100
#test_proteinsplit("dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@split="+str(split))

vsplit=20
test_proteinvsplit("variant=*@dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@vsplit="+str(vsplit))
