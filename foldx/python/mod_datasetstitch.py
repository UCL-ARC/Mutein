"""
RSA 11.4.22

Script to take a file in the format given to me by Michael Hall and then produce a list of proteins

We use the bioservices python library to access databases
https://pypi.org/project/bioservices/

"""
import os
import sys
import yaml
import pandas as pd
from os.path import exists
from shutil import copyfile

import _helper
import Paths
import Arguments
import FileDf


def run_pipeline(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    dataset_path = Paths.Paths(
        data_dir, install_dir, dataset=dataset,readonly=False
    )
    
    # load the list of the genes
    genes_fd = FileDf.FileDf(dataset_path.dataset_inputs + "genes_list.csv")
    genes_csv = genes_fd.openDataFrame()
        
    for g in range(len(genes_csv.index)):
        gn = genes_csv["gene"][g].lower()
        argsgn = args
        arglist = args[1]
        arglist += "@gene=" + gn
        argsgn[1] = arglist
        import mod_genestitch as ppi
        ppi.run_pipeline(argsgn)
        # now copy all available gene+results to dataset results
        gene_path = Paths.Paths(data_dir, install_dir, dataset=dataset,gene=gn,readonly=False)        
        filenameA = gene_path.outputs + "ddg_background.csv"
        filenameA_ds = dataset_path.outputs + f"{gn}_ddg_background.csv"
        filenameB = gene_path.outputs + "ddg_variants.csv"
        filenameB_ds = dataset_path.outputs + f"{gn}_ddg_variants.csv"
        if exists(filenameA):
            copyfile(filenameA, filenameA_ds)
        if exists(filenameB):
            copyfile(filenameB, filenameB_ds)

        
    #print("### COMPLETED dataset preparation ###")
    #print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
