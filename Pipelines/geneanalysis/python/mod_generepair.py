"""
-----------------------------
RSA 15/03/2022
-----------------------------
Adapted from
https://github.com/shorthouse-mrc/COVID_structure/blob/main/Foldx_repair6.py
-----------------------------

This file takes a pdb code (file must be located in the same directory in the format 1xyz.pdb)
It formats the pdb file into a paramater file suitable for foldx PositionScan
The output is in the same directory with the name
scanparams_1xyz.txt
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
import os
from shutil import copyfile


# import from the shared library in Mutein/Pipelines/shared/lib
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)
import Paths
import Arguments
import Config
import Foldx
import FileDf


def run_pipeline(args):

    ##### INPUTS #############################################
    # The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
    print("### FoldX repair job ###")
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")
    sys.path.append(install_dir)
    sys.path.append(install_dir + "/Pipelines")
    sys.path.append(install_dir + "/Pipelines/libs")
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    task = int(argus.arg("task", "none"))
    
    gene_path = Paths.Paths(        
        data_dir,
        install_dir + "Pipelines/geneanalysis",
        dataset=dataset,
        gene=gene,        
    )
    pdbtasks = gene_path.gene_outputs + "pdb_tasklist.csv"
    fio = FileDf.FileDf(pdbtasks)
    df = fio.openDataFrame()
    
    if task <= len(df.index):
        pdbcode = df["pdb"][task-1].lower()
        argsgn = args
        arglist = args[1]
        arglist += "@pdb=" + pdbcode
        argsgn[1] = arglist
        import Pipelines.geneanalysis.python.mod_pdbrepair as ppl
        print("Repairing pdb", pdbcode)
        ppl.run_pipeline(argsgn)
                
        
    else:
        print("Task beyond the data")

    print("### COMPLETED FoldX repair job ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
