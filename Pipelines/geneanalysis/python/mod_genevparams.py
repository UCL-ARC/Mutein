"""
-----------------------------
RSA 23/03/2022
-----------------------------
This file takes a pdb code (the data must be in inputs)
It formats the pdb file into a parameter file suitable for foldx PositionScan
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
import os
import pandas as pd
from shutil import copyfile
from os.path import exists

# import from the shared library in Mutein/Pipelines/shared/lib
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)
import Paths
import Arguments
import Config
import FileDf


##### INPUTS #############################################
# The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
# This version of the script does not do commbinations only the mutations
def run_pipeline(args):
    print("### FoldX make variant params job ###")
    print(args)
    ##############################################
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")
    sys.path.append(install_dir)
    sys.path.append(install_dir + "/Pipelines")
    sys.path.append(install_dir + "/Pipelines/libs")
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")

    gene_path = Paths.Paths(
        data_dir,
        install_dir + "Pipelines/geneanalysis",
        dataset=dataset,
        gene=gene,
    )
    pdbtasks = gene_path.gene_outputs + "pdb_tasklist.csv"
    fio = FileDf.FileDf(pdbtasks)
    df = fio.openDataFrame()

    all_params = []

    for t in range(len(df.index)):
        pdbcode = df["pdb"][t].lower()

        pdb_path = Paths.Paths(
            data_dir,
            install_dir + "Pipelines/geneanalysis",
            dataset=dataset,
            gene=gene,
            pdb=pdbcode,
        )
        argsgn = args
        arglist = args[1]
        arglist += "@pdb=" + pdbcode
        argsgn[1] = arglist
        import Pipelines.geneanalysis.python.mod_pdbvparams as ppl

        ppl.run_pipeline(argsgn)
        filename = pdb_path.pdb_thruputs + "params_variants.txt"
        if exists(
            filename
        ):  # TODO we might want to make it optional to fail if there is a missing file
            fdf = FileDf.FileDf(filename, sep=" ")
            csv = fdf.openDataFrame()
            all_params.append(csv)

    all_path = gene_path.gene_thruputs + "params_variants.txt"
    if len(all_params) > 0:
        all_df = pd.concat(all_params, axis=0)
        all_df.to_csv(all_path, index=False, sep=" ", header=True)
        print("MUTEIN SCRIPT ENDED")
        return True
    print("MUTEIN SCRIPT ENDED")
    return False


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
