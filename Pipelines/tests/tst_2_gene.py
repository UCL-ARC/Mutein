"""
RSA: 20/5/22
This CI test script runs from the gene level
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
def test_geneprep(inputs):
    inputs = addpath(inputs)
    import ga_2_0_geneprep as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_genesplit(inputs):
    inputs = addpath(inputs)
    import ga_2_genebackparams as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_genevsplit(inputs):
    inputs = addpath(inputs)
    import ga_2_genevarparams as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_generepair(inputs):
    inputs = addpath(inputs)
    import ga_2_generepair as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_genetasks(inputs):
    inputs = addpath(inputs)
    import ga04a_posscan as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_genevtasks(inputs):
    inputs = addpath(inputs)
    import ga04b_singlescan as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

def test_genedblagg(inputs):
    inputs = addpath(inputs)
    import ga_2_genestitch as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)
######################################################################
### INPUTS
dataset=""
gene="notch1"

repairs=1
split=10000
vsplit=2000

# whhich steps of the pipeline to run
prepareA = 0
prepareB = 0
repair = 0
tasks = 0
vtasks = 0
doubleagg = 1

if prepareA:
    # @@@@ - PREPARE A - @@@@
    test_geneprep("variant=*@dataset="+dataset+"@gene="+gene+"@split="+str(split)+"@vsplit="+str(vsplit))
if prepareB:
    # @@@@ - PREPARE B - @@@@
    test_genesplit("dataset="+dataset+"@gene="+gene+"@split="+str(split))
    test_genevsplit("variant=*@dataset="+dataset+"@gene="+gene+"@vsplit="+str(vsplit))
if repair:
    # @@@@ - REPAIR - @@@@
    test_generepair("dataset="+dataset+"@gene="+gene+"@repairs="+str(repairs)+"@task=2")
if tasks:
    # @@@@ - TASKS - @@@@
    test_genetasks("dataset="+dataset+"@gene="+gene+"@repairs="+str(repairs)+"@task=2636")
if vtasks:
    # @@@@ - Variant TASKS - @@@@
    test_genevtasks("dataset="+dataset+"@gene="+gene+"@repairs="+str(repairs)+"@task=1150")
if doubleagg:
    # @@@@ - Variant AGG - @@@@
    test_genedblagg("dataset="+dataset+"@gene="+gene+"@repairs="+str(repairs))
