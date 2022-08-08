"""
-----------------------------
RSA 15/03/2022
-----------------------------
This file takes a pdb code (file must be located in the same directory in the format 1xyz.pdb)
It formats the pdb file into a parameter file suitable for foldx PositionScan
The output is in the same directory with the name
scanparams_1xyz.txt
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
# import sys
import os
import pandas as pd
from os.path import exists
import sys

import _helper
import Paths
import Arguments
import Config
import FileDf

##### INPUTS #############################################
# The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
def run_pipeline(args):
    print("### FoldX make params job ###")
    print(args)
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    split = int(argus.arg("split", 0))

    gene_path = Paths.Paths(
        data_dir,
        install_dir,
        dataset=dataset,
        gene=gene,
    )
    pdbtasks = gene_path.gene_outputs + "pdb_tasklist.csv"
    if not exists(pdbtasks):
        #print("MUTEIN SCRIPT ENDED")
        return False

    fio = FileDf.FileDf(pdbtasks)
    df = fio.openDataFrame()

    all_params = []

    numpdbs = len(df.index)
    # There is a 100,000 limit on number of tasks, so if the pdbs make more than that, reduce number of tasks
    # This could be a problem fir gene AURKA which has many pdbs
    if split == 0:
        chunk = int(argus.arg("chunk"), 0)
        split = int(numpdbs / chunk)
    if numpdbs * split > 100000:
        split = int(100000 / numpdbs) - 1
        args[1] += "@split=" + str(split)

    for t in range(len(df.index)):
        pdbcode = df["pdb"][t].lower()
        if pdbcode[0] != "#":
            pdb_path = Paths.Paths(
                data_dir,
                install_dir,
                dataset=dataset,
                gene=gene,
                pdb=pdbcode,
            )
        

            argsgn = args
            arglist = args[1]
            arglist += "@pdb=" + pdbcode
            argsgn[1] = arglist
            import mod_pdbparams as ppl

            ppl.run_pipeline(argsgn)
            filename = pdb_path.pdb_inputs + "params_background.txt"
            if exists(
                filename
            ):  # TODO we might want to make it optional to fail if there is a missing file
                fdf = FileDf.FileDf(filename, sep=" ")
                csv = fdf.openDataFrame()
                all_params.append(csv)

    all_path = gene_path.gene_inputs + "params_background.txt"
    if len(all_params) > 0:
        all_df = pd.concat(all_params, axis=0)
        all_df.to_csv(all_path, index=False, sep=" ", header=True)
        #print("MUTEIN SCRIPT ENDED")
        return True
    #print("MUTEIN SCRIPT ENDED")
    return False


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
