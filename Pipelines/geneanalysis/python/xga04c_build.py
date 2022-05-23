"""
-----------------------------
RSA 15/03/22
-----------------------------
Adapted from:
https://github.com/shorthouse-mrc/COVID_structure/blob/main/Alpha_variant/Foldx_variantcombinations_calculation.py
-----------------------------

This performs a mutation on a given list of amino acid positions on the structure
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
import os
import pandas as pd
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


def run_pipeline06(args):
    print("### FoldX build job ###")
    print(args)
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")
    sys.path.append(install_dir)
    sys.path.append(install_dir + "/Pipelines")
    sys.path.append(install_dir + "/Pipelines/libs")
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    pdbcode = argus.arg("pdb").lower()
    pdb_path = Paths.Paths(
        "pdb",
        data_dir,
        install_dir + "Pipelines/geneanalysis",
        dataset=dataset,
        gene=gene,
        pdb=pdbcode,
    )
    # pdb_config = Config.Config(pdb_path.pdb_inputs + "/config.yml")
    # argus.addConfig(pdb_config.params)
    task = argus.arg("task", "none")
    mutation_string = argus.arg("mutation", "none")

    ############################################
    # set up the files and directories
    pdb = pdbcode + "_rep" + str(argus.arg("repairs"))
    pdbfile = pdbcode + "_rep" + str(argus.arg("repairs")) + ".pdb"
    mutations = []

    if mutation_string == "none":
        filename = pdb_path.pdb_thruputs +"params_variants.txt"
        print("open", filename)
        with open(filename) as fr:
            paramscontent = fr.readlines()
            if task == "all":
                for row in paramscontent:
                    row = row.strip()
                    print(row)
                    rowvals = row.split(" ")
                    mutation = rowvals[2]
                    row = rowvals[3]
                    mutations.append([mutation, row])
            else:
                row = paramscontent[int(task) - 1].strip()
                rowvals = row.split(" ")
                mutation = rowvals[2]
                row = rowvals[3]
                mutations.append([mutation, row])
    else:
        # we have specified a mutation and row from the file
        mutations.append([mutation_string, 0])

    for mut, row in mutations:
        print(mut, row)

        row_path = (
            pdb_path.pdb_thruputs
            + str(argus.arg("split"))
            + "_"
            + str(row)
            + "_build"
            + "/"
        )
        print("### ... change directory", row_path)
        argus.params["thisrow"] = row
        argus.params["thismut"] = mut
        pdb_path.goto_job_dir(row_path, args, argus.params, "_inputs06")
        print(
            "### foldx06: ... copying file",
            pdb_path.pdb_thruputs + pdbfile,
            row_path + pdbfile,
        )
        copyfile(pdb_path.pdb_thruputs + pdbfile, row_path + pdbfile)

        fx_runner = Foldx.Foldx(argus.arg("foldxe"))
        #fx_runner.runBuild(pdbfile, mut, 15)
        #Dif_af-p46531-f1-model_v2_rep1.fxout
        ddg_file = row_path + "Dif_" + pdb + ".fxout"
        df_file = row_path + "pdbfile_build_DDG.csv"
        
        filename = pdb_path.pdb_inputs + "coverage.csv"
        fdfp = FileDf.FileDf(filename)
        cov_df = fdfp.openDataFrame()

        #path, pdbfile, pdb_mut, gene_mut, coverage, outfile_path
        #fx_runner.createBuildCsv(pdbcode, mut, ddg_file, df_file, task)
        fx_runner.createBuildCsv(row_path,pdbcode, mut, ddg_file, df_file, task)


if __name__ == "__main__":
    import sys

    globals()["run_pipeline06"](sys.argv)
