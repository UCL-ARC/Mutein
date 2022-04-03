
import sys,os

def addpath():
    above_path = os.path.dirname(os.path.realpath(__file__))[:-5]    
    sys.path.append(above_path)

def test_foldx_pipeline_HIS1():
    addpath()
    import foldx03_posscan as ppl
    args = ['','pdb=6vxx','name=6vxx_tst','mutation=VA47a,VA47a','row=005', 'repairs=5']
    ppl.run_pipeline03(args)

def test_foldx_pipeline_HIS2():
    addpath()
    import foldx03_posscan as ppl
    args = ['','pdb=6vxx','name=6vxx_tst','mutation=VA47a,VA47a','row=010', 'repairs=10']
    ppl.run_pipeline03(args)
#########################################################################
test_foldx_pipeline_HIS1()
test_foldx_pipeline_HIS2()