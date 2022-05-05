"""
RSA: 2/4/22
This CI test script runs each script in the pipleline with a very small dataset
It enables debugging of the scripts as if run from a batch

"""
import sys, os
# import from the shared library in Mutein/Pipelines/shared/lib
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/shared/libs"
sys.path.append(retpath)
retpath = "/".join(dirs) + "/foldx/scripts"
sys.path.append(retpath)
retpath = "/".join(dirs) + "/geneprot/scripts"
sys.path.append(retpath)
retpath = "/".join(dirs) + "/genestitch"
sys.path.append(retpath)
import Paths

def addpath():        
    print(sys.path)

######################################################################
def test_geneprot_pipeline(inputs):
    addpath()
    import pipeline02_genestoproteins as p02
    args = ["", inputs]
    p02.run_pipeline(args)

def test_genestitch_ppl(inputs):
    addpath()
    import pipeline01_genestitch as p01
    args = ["", inputs]
    p01.run_pipeline(args)

########################################################################
def test_foldx_pipeline_repair(inputs):
    addpath()
    import stitch01_genestitch as p01
    args = ["", inputs]
    p01.run_pipeline01(args)

def test_foldx_pipeline_params(inputs):
    addpath()
    import foldx02_makeparams as p02
    args = ["", inputs]
    p02.run_pipeline02(args)


def test_foldx_pipeline_posscan(inputs):
    addpath()
    import foldx03_posscan as p03
    args = ["", inputs]
    p03.run_pipeline03(args)

def test_foldx_pipeline_agg(inputs):
    addpath()
    import foldx04_aggddg as p04
    args = ["", inputs]
    p04.run_pipeline04(args)

# variants
def test_foldx_pipeline_vparams(inputs):
    addpath()
    import foldx05_singlesparams as p05
    args = ["", inputs]
    p05.run_pipeline05(args)

def test_foldx_pipeline_build(inputs):
    addpath()
    import foldx06_build as p06
    args = ["", inputs]
    p06.run_pipeline06(args)

def test_foldx_pipeline_singlescan(inputs):
    addpath()
    import foldx06_singlescan as p06
    args = ["", inputs]
    p06.run_pipeline03(args)

def test_foldx_pipeline_singlesagg(inputs):
    addpath()
    import foldx07_singlesagg as p07
    args = ["", inputs]
    p07.run_pipeline04(args)

######################################################################

## GENE_PROT
#test_geneprot_pipeline("dataset=notch@")

## FOLDX
#test_foldx_pipeline_repair("repairs=2@dataset=notch@gene=NOTCH1@pdb=3v79@chain=A")

#test_foldx_pipeline_params("repairs=2@split=200@dataset=notch@gene=NOTCH1@pdb=3v79@chain=RK@task=1")
#test_foldx_pipeline_posscan("repairs=2@split=200@dataset=notch@gene=NOTCH1@pdb=3v79@task=2")
#test_foldx_pipeline_agg("repairs=2@split=200@dataset=notch@gene=NOTCH1@pdb=3v79")

#test_foldx_pipeline_vparams("repairs=2@split=20@dataset=notch@gene=NOTCH1@pdb=3v79@task=1@variant=*")
###test_foldx_pipeline_build("repairs=2@split=20@dataset=notch@gene=NOTCH1@pdb=3v79@task=1")
#test_foldx_pipeline_singlescan("repairs=2@split=20@dataset=notch@gene=NOTCH1@pdb=3v79@@task=1")
#test_foldx_pipeline_singlescan("repairs=2@split=20@dataset=notch@gene=NOTCH1@pdb=3v79@@task=5")
#test_foldx_pipeline_singlesagg("repairs=2@split=20@dataset=notch@gene=NOTCH1@pdb=3v79")

test_genestitch_ppl("dataset=notch@gene=NOTCH1")
