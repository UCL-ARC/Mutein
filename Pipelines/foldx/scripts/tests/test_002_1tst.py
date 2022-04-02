
'''
RSA: 2/4/22
This CI test script runs the entire pipeline with a very small dataset
'''
import sys,os

def addpath():
    above_path = os.path.dirname(os.path.realpath(__file__))[:-5]    
    sys.path.append(above_path)

def test_foldx_pipeline01_1tst():
    addpath()
    import foldx01_repair_mod as p01
    args = ['','user=CI','pdb=1tst']    
    p01.run_pipeline01(args)
    
def test_foldx_pipeline02_1tst():
    addpath()
    import foldx02_makeparams_mod as p02
    args = ['','user=CI','pdb=1tst']    
    p02.run_pipeline02(args)

def test_foldx_pipeline03_1tst():
    addpath()
    import foldx03_posscan_mod as p03
    args = ['','user=CI','pdb=1tst','mutation=.']    
    p03.run_pipeline03(args)

def test_foldx_pipeline04_1tst():
    addpath()
    import foldx04_aggddg_mod as p04
    args = ['','user=CI','pdb=1tst']        
    p04.run_pipeline04(args)

test_foldx_pipeline01_1tst()
test_foldx_pipeline02_1tst()
test_foldx_pipeline03_1tst()
test_foldx_pipeline04_1tst()
