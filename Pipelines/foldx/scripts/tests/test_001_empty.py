
import sys,os

def addpath():
    above_path = os.path.dirname(os.path.realpath(__file__))[:-5]    
    sys.path.append(above_path)

def test_foldx_pipeline00_emptyA():
    addpath()
    import foldx00_pipeline as p00
    args = ['','user=empty_hpc','pdb=1tst','jobs=1234','id=1:time=6:00:00']    
    p00.run_pipeline00(args)

def test_foldx_pipeline00_emptyB():
    addpath()
    import foldx00_pipeline as p00
    args = ['','user=empty_hpc','pdb=6vxx','jobs=1234','id=1@time=0:05:00','id=3@time=2:00:00','id=3@array=200']
    p00.run_pipeline00(args)

#test_foldx_pipeline00_emptyA()
test_foldx_pipeline00_emptyB()
