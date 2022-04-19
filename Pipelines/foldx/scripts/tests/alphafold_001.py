"""
RSA 19/4/22
This checks the inputs and simple single run work on an alpha fold structure
"""
import sys, os


def addpath():
    above_path = os.path.dirname(os.path.realpath(__file__))[:-5]
    sys.path.append(above_path)


def test_foldx_pipeline_1():
    addpath()
    import foldx01_repair as ppl
    args = [
        "",
        "pdb=af-p46531-f1-model_v2",                        
        "repairs=1"
    ]
    ppl.run_pipeline01(args)

def test_foldx_pipeline_2():
    addpath()
    import foldx02_makeparams as ppl
    args = [
        "",
        "pdb=af-p46531-f1-model_v2",                        
        "repairs=1",
        "split=500"
    ]
    ppl.run_pipeline02(args)

def test_foldx_pipeline_3():
    addpath()
    import foldx03_posscan as ppl

    args = [
        "",
        "pdb=af-p46531-f1-model_v2",                
        "row=1",
        "repairs=1",
    ]
    ppl.run_pipeline03(args)


#########################################################################
#test_foldx_pipeline_1()
#test_foldx_pipeline_2()
test_foldx_pipeline_3()

