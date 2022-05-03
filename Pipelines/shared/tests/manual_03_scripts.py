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
import Paths

def addpath():        
    print(sys.path)

def test_foldx_pipeline_repair(inputs):
    addpath()
    import foldx01_repair as p01
    args = ["", inputs]
    p01.run_pipeline01(args)

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

#test_foldx_pipeline_repair("repairs=2@dataset=notch@gene=NOTCH1@pdb=1toz@chain=A")

#test_foldx_pipeline_posscan("repairs=5@split=100@dataset=notch@gene=NOTCH1@pdb=1toz@chain=A@task=1")
#test_foldx_pipeline_posscan("repairs=5@split=100@dataset=notch@gene=NOTCH1@pdb=1toz@chain=A@task=2")
test_foldx_pipeline_agg("repairs=5@split=100@dataset=notch@gene=NOTCH1@pdb=1toz@chain=A")

#test_foldx_pipeline_build("repairs=5@split=20@dataset=notch@gene=NOTCH1@pdb=1toz@chain=A@task=1")
#test_foldx_pipeline_singlescan("repairs=5@split=20@dataset=notch@gene=NOTCH1@pdb=1toz@chain=A@task=1")
#test_foldx_pipeline_singlescan("repairs=5@split=20@dataset=notch@gene=NOTCH1@pdb=1toz@chain=A@task=2")
test_foldx_pipeline_singlesagg("repairs=5@split=20@dataset=notch@gene=NOTCH1@pdb=1toz@chain=A")
