"""
RSA: 2/4/22
This CI test script runs each script in the pipleline with a very small dataset

"""
import sys, os


def addpath():
    above_path = os.path.dirname(os.path.realpath(__file__))[:-5]
    sys.path.append(above_path)


def test_foldx_pipeline01_1tst():
    addpath()
    import foldx01_repair as p01
    args = ["", "user=CI", "pdb=1tst","repairs=3"]
    p01.run_pipeline01(args)

def test_foldx_pipeline02_1tst():
    addpath()
    import foldx02_makeparams as p02
    args = ["", "user=CI", "pdb=1tst","repairs=3"]
    p02.run_pipeline02(args)


def test_foldx_pipeline03_1tst():
    addpath()
    import foldx03_posscan as p03
    args = ["", "user=CI", "pdb=1tst", "row=.","repairs=3"]
    p03.run_pipeline03(args)

def test_foldx_pipeline04_1tst():
    addpath()
    import foldx04_aggddg as p04
    args = ["", "user=CI", "pdb=1tst","repairs=3"]
    p04.run_pipeline04(args)

def test_foldx_pipeline05_1tst():
    addpath()
    import Pipelines.foldx.scripts.foldx05_buildparams as p05
    args = ["", "user=CI", "pdb=1tst","repairs=3"]
    p05.run_pipeline05(args)


#test_foldx_pipeline01_1tst()
#test_foldx_pipeline02_1tst()
#test_foldx_pipeline03_1tst()
#test_foldx_pipeline04_1tst()
test_foldx_pipeline05_1tst()
