
import sys,os

def addpath():
    above_path = os.path.dirname(os.path.realpath(__file__))[:-5]    
    sys.path.append(above_path)

def test_foldx_pipeline_tsta():
    addpath()
    import foldx03_posscan as ppl
    args = ['','pdb=1tst','name=1tst_2','mutation=TA29a']
    ppl.run_pipeline03(args)

def test_foldx_pipeline_tstb():
    addpath()
    import foldx03_posscan as ppl
    args = ['','pdb=1tst','name=1tst_2','mutation=TA29a','row=row1']
    ppl.run_pipeline03(args)

def test_foldx_pipeline_tstc():
    addpath()
    import foldx03_posscan as ppl
    args = ['','pdb=1tst','name=1tst_2','row=2']
    ppl.run_pipeline03(args)

def test_foldx_pipeline_tstall():
    addpath()
    import foldx03_posscan as ppl
    args = ['','pdb=1tst','name=1tst_2','row=.']
    ppl.run_pipeline03(args)


#########################################################################
#test_foldx_pipeline_tsta()
#test_foldx_pipeline_tstb()
#test_foldx_pipeline_tstc()
test_foldx_pipeline_tstall()