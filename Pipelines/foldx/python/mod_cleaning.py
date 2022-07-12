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

def deletePath(thispath,remPath=False):
    print("Cleaning up",thispath)
    fs = []
    ps = []
    for (path, dirnames, filenames) in os.walk(thispath):
        fs.extend(filenames)
        ps.extend(dirnames)    
    for f in fs:                
        onepath = thispath + f
        #print(onepath)
        if os.path.exists(onepath):
            os.remove(onepath)                
    for p in ps:
        if p.upper() != "RESULTS":
            onepath = thispath + p
            #print(onepath)
            shutil.rmtree(onepath,ignore_errors=True)
    if remPath:
        shutil.rmtree(thispath,ignore_errors=True)
            
###############################################################################
def cleanGene(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")    
    gene_path = Paths.Paths(data_dir, install_dir, dataset=dataset,gene=gene,readonly=False)    
    mypath = gene_path.gene_inputs    
    deletePath(mypath)

def cleanGeneThruputs(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")    
    gene_path = Paths.Paths(data_dir, install_dir, dataset=dataset,gene=gene,readonly=False)    
    mypath = gene_path.gene_thruputs
    deletePath(mypath)

def cleanPdbThruputs(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")    
    pdb = argus.arg("pdb")    
    pdb_path = Paths.Paths(data_dir, install_dir, dataset=dataset,gene=gene,pdb=pdb,readonly=False)    
    mypath = pdb_path.pdb_thruputs
    deletePath(mypath)

def cleanPdbThruBack(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")    
    pdb = argus.arg("pdb")    
    task = argus.arg("task")
    pdb_path = Paths.Paths(data_dir, install_dir, dataset=dataset,gene=gene,pdb=pdb,readonly=False)    
    mypath = pdb_path.pdb_thruputs + "back_" + str(task)
    deletePath(mypath,remPath=True)

def cleanPdbThruVar(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")    
    pdb = argus.arg("pdb")    
    task = argus.arg("task")
    pdb_path = Paths.Paths(data_dir, install_dir, dataset=dataset,gene=gene,pdb=pdb,readonly=False)    
    mypath = pdb_path.pdb_thruputs + "var_" + str(task)
    deletePath(mypath,remPath=True)

    
    


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
