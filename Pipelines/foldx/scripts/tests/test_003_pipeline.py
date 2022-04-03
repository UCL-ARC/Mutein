
import sys,os

def addpath():
    above_path = os.path.dirname(os.path.realpath(__file__))[:-5]    
    sys.path.append(above_path)

def test_foldx_pipeline00_1tstA():
    addpath()
    import foldx00_pipeline as p00
    args = ['','user=CI','jobs=1234','pdb=1tst']    
    p00.run_pipeline00(args)

def test_foldx_pipeline00_1tstB():
    addpath()
    import foldx00_pipeline as p00
    args = ['','user=CI','jobs=1234','pdb=1tst','repairs=15']    
    p00.run_pipeline00(args)

#test_foldx_pipeline00_1tstA()
test_foldx_pipeline00_1tstB()
