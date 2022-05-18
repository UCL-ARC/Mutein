"""
RSA: 2/4/22
This CI test script runs each script in the pipleline with a very small dataset
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
def test_geneprot_pipeline(inputs):
    inputs = addpath(inputs)
    import ga01_genestoproteins as p02
    args = ["", inputs]
    p02.run_pipeline(args)

def test_genestitch_ppl(inputs):
    inputs = addpath(inputs)
    import ga06_genestitch as p01
    args = ["", inputs]
    p01.run_pipeline(args)

########################################################################
def test_foldx_pipeline_repair(inputs):
    inputs = addpath(inputs)
    import ga02_repair as p01
    args = ["", inputs]
    p01.run_pipeline01(args)

def test_foldx_pipeline_params(inputs):
    inputs = addpath(inputs)
    import ga03a_makeparams as p02
    args = ["", inputs]
    p02.run_pipeline02(args)


def test_foldx_pipeline_posscan(inputs):
    inputs = addpath(inputs)
    import ga04a_posscan as p03
    args = ["", inputs]
    p03.run_pipeline03(args)

def test_foldx_pipeline_agg(inputs):
    inputs = addpath(inputs)
    import ga05a_aggddg as p04
    args = ["", inputs]
    p04.run_pipeline04(args)

# variants
def test_foldx_pipeline_vparams(inputs):
    inputs = addpath(inputs)
    import ga03b_singlesparams as p05
    args = ["", inputs]
    p05.run_pipeline05(args)

def test_foldx_pipeline_build(inputs):
    inputs = addpath(inputs)
    import ga04c_build as p06
    args = ["", inputs]
    p06.run_pipeline06(args)

def test_foldx_pipeline_singlescan(inputs):
    inputs = addpath(inputs)
    import ga04b_singlescan as p06
    args = ["", inputs]
    p06.run_pipeline03(args)

def test_foldx_pipeline_singlesagg(inputs):
    inputs = addpath(inputs)
    import ga05b_singlesagg as p07
    args = ["", inputs]
    p07.run_pipeline04(args)

######################################################################
### INPUTS
dataset="any"
gene = "any"
repairs = 1
split=100
vsplit=20
pdbtask = 2

test_geneprot_pipeline("dataset="+dataset+"@")
for pdbcode in ["x"]:        
    print(pdbcode)
    pdb_add = ""
    if pdbcode != "x":
        pdb_add = "@pdb="+pdbcode
    #test_foldx_pipeline_repair("repairs="+str(repairs) + "@dataset="+dataset+"@gene="+gene+"@task="+str(pdbtask))
    #############test_foldx_pipeline_background
    #test_foldx_pipeline_params("repairs="+str(repairs)+"@split="+str(split)+"@dataset="+dataset+"@gene="+gene)
    #test_foldx_pipeline_posscan("repairs="+str(repairs)+"@dataset="+dataset+"@gene="+gene+"@task=1001")
    #test_foldx_pipeline_posscan("repairs="+str(repairs)+"@split="+str(split)+"@dataset="+dataset+"@gene="+gene+"@pdb="+pdbcode+"@task=2")
    #test_foldx_pipeline_posscan("repairs="+str(repairs)+"@split="+str(split)+"@dataset="+dataset+"@gene="+gene+"@pdb="+pdbcode+"@task=3")
    #test_foldx_pipeline_agg("repairs="+str(repairs)+"@dataset="+dataset+"@gene="+gene+"@pdb="+pdbcode)
    #############test_foldx_pipeline_variants
    #test_foldx_pipeline_vparams("repairs="+str(repairs)+"@split="+str(vsplit)+"@dataset="+dataset+"@gene="+gene+"@variant=*")
    #test_foldx_pipeline_singlescan("repairs="+str(repairs)+"@dataset="+dataset+"@gene="+gene+"@pdb="+pdbcode+"@task=154")
    #test_foldx_pipeline_singlescan("repairs="+str(repairs)+"@split="+str(vsplit)+"@dataset="+dataset+"@gene="+gene+"@pdb="+pdbcode+"@task=2")
    #test_foldx_pipeline_singlescan("repairs="+str(repairs)+"@split="+str(vsplit)+"@dataset="+dataset+"@gene="+gene+"@pdb="+pdbcode+"@task=3")
    #test_foldx_pipeline_singlesagg("repairs="+str(repairs)+"@dataset="+dataset+"@gene="+gene+pdb_add)
## GENE_STITCH
#test_genestitch_ppl("dataset=notch@gene=NOTCH1")
