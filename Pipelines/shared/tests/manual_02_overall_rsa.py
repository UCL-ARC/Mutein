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
    


def test_pipeline(method):
    addpath()
    import overall_rsa as pipe
    args = [
        "",
        "batch_tst01.yml",                        
        method
    ]
    pipe.overall_rsa(args)


######################################
test_pipeline("qsub")
#test_pipeline("py")
#test_pipeline("sh")
