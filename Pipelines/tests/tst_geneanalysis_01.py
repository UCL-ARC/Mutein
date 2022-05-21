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
        
def test_pipeline(method, batch_file,dataset,gene,pdb):
    addpath()
    # There are 7 arguments
    #install_dir = args[1]   # 1) the executable installation directory, the root directory of the peipeline    
    #working_dir = args[2]   # 1) the working dir, the root that the data output and input lives in
    #yaml_file = args[3]     # 2) a yaml file path with the batch definition
    #py_or_sh = args[4]      # 3) qsub or py or sh for python or hpc batch or just sh        
    #dataset = args[5]  # 3) dataset    
    #gene = args[6]  # 4) gene    
    #pdb = args[7]  # 5) pdb

    import Pipelines.libs.pipeline_qsubber as pipe
    args = [
        "",
        "/home/rachel/UCL/github/Mutein/",
        "/home/rachel/UCL/github/MuteinData/",
        "/home/rachel/UCL/github/Mutein/Pipelines/geneanalysis/config/" + batch_file,                           
        method,
        dataset,
        gene,
        pdb
    ]
    pipe.pipeline_qsubber(args)


######################################
## Tests for the geneprot
#test_pipeline("py","batch_geneprot.yml","notch","NOTCH1")

## Tests for the pdb
#test_pipeline("py","batch_pdb.yml","notch","NOTCH1","1toz","A")
#test_pipeline("qsub","batch_tst02.yml")
#test_pipeline("qsub_tst","batch_tst02.yml")
#test_pipeline("py","batch_tst02.yml")
#test_pipeline("qsub_tst","batch_notch_alpha.yml")
#test_pipeline("sh","batch_tst02.yml")

## Tests for the gene stitch
#test_pipeline("qsub_tst","batch_gene_tasks.yml","notch","NOTCH1","")
test_pipeline("py","batch_dataset_tasks.yml","notch","","")
#test_pipeline("py","batch_pdb_agg.yml","","","1pb5")


