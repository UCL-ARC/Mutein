
import sys,os

def addpath():
    above_path = os.path.dirname(os.path.realpath(__file__))[:-5]    
    sys.path.append(above_path)

def test_foldx_pipeline00_empty():
    addpath()
    import foldx00_pipeline_mod as p00
    args = ['','user=empty']    
    p00.run_pipeline00(args)

def test_foldx_pipeline01_test():
    addpath()
    import foldx01_repair_mod as p01
    args = ['','user=CI','pdb=Test','jobs=1']    
    p01.run_pipeline01(args)
    
#test_foldx_pipeline01_test()