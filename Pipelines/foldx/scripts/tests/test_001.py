
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
    
def test_foldx_pipeline02_test():
    addpath()
    import foldx02_makeparams_mod as p02
    args = ['','user=CI','pdb=Test']    
    p02.run_pipeline02(args)

def test_foldx_pipeline02_covid():
    addpath()
    import foldx02_makeparams_mod as p02
    args = ['','user=CI','pdb=6vxx','split=25','name=covid1']    
    p02.run_pipeline02(args)
    
def test_foldx_pipeline03_test():
    addpath()
    import foldx03_posscan_mod as p03
    args = ['','user=CI','pdb=Test','mutation=AA27a,YA28a,TA29a','row=1']    
    p03.run_pipeline03(args)

def test_foldx_pipeline03_covid():
    addpath()
    import foldx03_posscan_mod as p03
    args = ['','user=CI','pdb=6vxx','mutation=AA27a','name=covid1']        
    p03.run_pipeline03(args)

def test_foldx_pipeline04_test():
    addpath()
    import foldx04_aggddg_mod as p04
    args = ['','user=CI','pdb=Test']        
    p04.run_pipeline04(args)

def test_foldx_pipeline05_test():
    addpath()
    import foldx05_vparams_mod as p05
    args = ['','user=CI','pdb=Test']        
    p05.run_pipeline05(args)

#test_foldx_pipeline00_empty()
#test_foldx_pipeline01_test()
#test_foldx_pipeline02_test()
#test_foldx_pipeline02_covid()
#test_foldx_pipeline03_test()
#test_foldx_pipeline03_covid()
#test_foldx_pipeline04_test()
#test_foldx_pipeline05_test()