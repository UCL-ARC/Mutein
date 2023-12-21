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
import sys

import _helper
import Paths
import Arguments
import Config
import Foldx
import FileDf


def run_pipeline(args):

    ##### INPUTS #############################################
    # The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
    print("### FoldX gene stitch job ###")
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    pdbcode = argus.arg("pdb", "").lower()
    gene_stitch = pdbcode == ""

    gene_path = Paths.Paths(
        data_dir,
        install_dir,
        dataset=dataset,
        gene=gene,
        pdb=pdbcode,
    )

    pdb_list = []

    if pdbcode != "":
        pdb_list.append(pdbcode)
    else:
        pdbtasks = gene_path.gene_outputs + "pdb_tasklist.csv"
        print("Open file:", pdbtasks)
        fio = FileDf.FileDf(pdbtasks)
        df = fio.openDataFrame()

        for t in range(len(df.index)):
            pdbcode = df["pdb"][t].lower()
            pdb_list.append(pdbcode)

    for pdbcode in pdb_list:
        if pdbcode[0] != "#": #we can comment out failed pdbs
            argsgn = []
            for arg in args:
                argsgn.append(arg)
            arglist = args[1]
            arglist += "@pdb=" + pdbcode
            argsgn[1] = arglist
            import mod_pdbagg as ppla

            print("Aggregating background pdb", pdbcode)
            ppla.run_pipeline(argsgn)
            import mod_pdbvagg as pplb

            print("Aggregating variants pdb", pdbcode)
            pplb.run_pipeline(argsgn)

    if gene_stitch:
        print("Gene stitching...")
        import mod_onegenestitch as pplc

        print("Stitching gene", gene)
        pplc.run_pipeline(args)
    else:
        print("No gene stitch, pdb only")

    print("### COMPLETED FoldX stitch job ###")
    #print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
