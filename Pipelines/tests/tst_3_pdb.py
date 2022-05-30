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
dataset="shearwater"

#gene="CR2"
#pdb="2gsx"
#pdb="1ly2"

gene="notch1"
pdb="1pb5"

repairs=8
split=20
vsplit=50

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
runs = "j"
runPPL("runs="+runs+"@dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@repairs="+str(repairs)+"@task=2"+"@chunk="+str(split)+"@variant=*@vchunk="+str(vsplit))


"""
# whhich steps of the pipeline to run
repair = 1
prepareA = 0
prepareB = 0
tasks = 0
vtasks = 0
agg = 0
vagg = 0
doubleagg = 0

if repair:
    # @@@@ - REPAIR - @@@@
    test_proteinrepair("dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@repairs="+str(repairs))
if prepareA:
    # @@@@ - PREPARE A - @@@@
    test_proteinprep("variant=*@dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@split="+str(split)+"@vsplit="+str(vsplit))
if prepareB:
    # @@@@ - PREPARE B - @@@@
    test_proteinsplit("dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@split="+str(split))
    test_proteinvsplit("variant=*@dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@vsplit="+str(vsplit))
if tasks:
    # @@@@ - TASKS - @@@@
    test_proteintasks("dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@repairs="+str(repairs)+"@task=1")
if vtasks:
    # @@@@ - Variant TASKS - @@@@
    test_proteinvtasks("dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@repairs="+str(repairs)+"@task=1")
if agg:
    # @@@@ - AGG - @@@@
    test_proteinagg("dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@repairs="+str(repairs))
if vagg:
    # @@@@ - Variant AGG - @@@@
    test_proteinvagg("dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@repairs="+str(repairs))
if doubleagg:
    # @@@@ - Variant AGG - @@@@
    test_proteindblagg("dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@repairs="+str(repairs))
"""