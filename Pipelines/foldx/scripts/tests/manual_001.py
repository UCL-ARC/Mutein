
import sys,os

def addpath():
    above_path = os.path.dirname(os.path.realpath(__file__))[:-5]    
    sys.path.append(above_path)

def test_foldx_pipeline_tst():
    addpath()
    import foldx03_posscan as ppl
    args = ['','pdb=1tst','name=1tst_2','mutation=TA29a']
    ppl.run_pipeline03(args)

#########################################################################
test_foldx_pipeline_tst()