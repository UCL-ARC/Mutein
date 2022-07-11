"""
RSA 11.4.22

Script to take a file in the format given to me by Michael Hall and then produce a list of proteins

We use the bioservices python library to access databases
https://pypi.org/project/bioservices/

"""
import os
import sys
import shutil

import _helper
import Paths
import Arguments
import FileDf


def run_pipeline(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")    
    
    gene_path = Paths.Paths(
        data_dir, install_dir, dataset=dataset,gene=gene,readonly=False
    )

    # walk through and get the directories and paths
    mypath = gene_path.gene_inputs
    fs = []
    ps = []
    for (path, dirnames, filenames) in os.walk(mypath):
        fs.extend(filenames)
        ps.extend(dirnames)
    
    for f in fs:                
        onepath = mypath + f
        #print(onepath)
        if os.path.exists(onepath):
            os.remove(onepath)
        
        
    for p in ps:
        if p.upper() != "RESULTS":
            onepath = mypath + p
            #print(onepath)
            shutil.rmtree(onepath,ignore_errors=True)
            
    
    


##########################################################################################
if __name__ == "__main__":
    import sys
    globals()["run_pipeline"](sys.argv)
