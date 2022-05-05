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
    


def test_pipeline(method, batch_file,dataset,gene,pdb,chain):
    addpath()
    import overall_rsa as pipe
    args = [
        "",
        batch_file,                           
        method,
        dataset,
        gene,
        pdb,
        chain
    ]
    pipe.overall_rsa(args)


######################################
## Tests for the geneprot
test_pipeline("qsub_tst","batch_genestitch.yml","notch","NOTCH1","","")

## Tests for the pdb
#test_pipeline("py","batch_pdb.yml","notch","NOTCH1","1toz","A")
#test_pipeline("qsub","batch_tst02.yml")
#test_pipeline("qsub_tst","batch_tst02.yml")
#test_pipeline("py","batch_tst02.yml")
#test_pipeline("qsub_tst","batch_notch_alpha.yml")
#test_pipeline("sh","batch_tst02.yml")

