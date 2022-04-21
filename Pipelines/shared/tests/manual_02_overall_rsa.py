# test_dummy.py
"""
RSA 21/4/2022
---------------------------
Any function thatstarts test_ in a file that startes test_ will be run automatically following a pull request.
This tests the outer pipeline
---------------------------

"""
import os
import sys

def addpath():        
    dirs = (os.path.dirname(os.path.realpath(__file__))).split("/")[:-2]
    newpath = "/".join(dirs)        
    print("Adding sys path=", newpath)
    sys.path.append(newpath)
    


def test_pipeline_sh():
    addpath()
    import overall_rsa as pipe
    args = [
        "",
        "batch_tst01.yml",                        
        "sh"
    ]
    pipe.overall_rsa(args)

def test_pipeline_py():
    addpath()
    import overall_rsa as pipe
    args = [
        "",
        "batch_tst01.yml",                        
        "py"
    ]
    pipe.overall_rsa(args)

######################################
test_pipeline_sh()
test_pipeline_py()