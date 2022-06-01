"""
-----------------------------
RSA 15/03/22
-----------------------------
Adapted from:
https://github.com/shorthouse-mrc/COVID_structure/blob/main/Foldx_positionscanall.py
-----------------------------

This performs a mutation on a given list of amino acid positions on the structure
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
# import sys
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
import Foldx
import FileDf


##### INPUTS #############################################
# The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)


def run_pipeline(args):
    print("### FoldX position scan job ###")
    print(args)
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")
    sys.path.append(install_dir)
    sys.path.append(install_dir + "/Pipelines")
    sys.path.append(install_dir + "/Pipelines/libs")
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    task = int(argus.arg("task", "0"))
    pdb = argus.arg("pdb", "")
    missing = argus.arg("missing", "N")

    if task == 0:
        print("No task entered")
    if task > 0:

        gene_path = Paths.Paths(
            data_dir,
            install_dir + "Pipelines/geneanalysis",
            dataset=dataset,
            gene=gene,
            pdb=pdb,
        )
        all_tasks = gene_path.thruputs + "params_background.txt"
        if missing.upper() == "Y":
            all_tasks = gene_path.thruputs + "params_background_incomplete.txt"
        print("Opening file", all_tasks)
        fio = FileDf.FileDf(all_tasks, sep=" ", header=True)
        df = fio.openDataFrame()

        if task <= len(df.index):
            pdbcode = df["pdb"][task - 1].lower()
            pdb_path = Paths.Paths(
                data_dir,
                install_dir + "Pipelines/geneanalysis",
                dataset=dataset,
                gene=gene,
                pdb=pdbcode,
            )
            # pdb_config = Config.Config(pdb_path.pdb_inputs + "/config.yml")
            # argus.addConfig(pdb_config.params)
            mutation_string = argus.arg("mutation", "none")
            ############################################
            pdbfile = pdbcode + "_rep" + str(argus.arg("repairs")) + ".pdb"
            pdbname = pdbcode + "_rep" + str(argus.arg("repairs"))
            mutations = []
            # task=all means all, task=1:n means an explicit row, row=-1 means the mutation string has been passd in explicitly
            if mutation_string == "none":
                if task == "all":
                    for i in range(len(df.index)):
                        mutation = df["mutation"][i]
                        row = df["row"][i]
                        mutations.append([mutation, row])
                else:
                    if int(task) <= len(df.index):
                        mutation = df["mutation"][int(task) - 1]
                        row = df["row"][int(task) - 1]
                        mutations.append([mutation, row])
            else:
                # we have specified a mutation and row from the file
                mutations.append([mutation_string, 0])

            work_path = pdb_path.pdb_thruputs + "agg/"
            pdb_path.goto_job_dir(work_path, args, argus.params, "_inputs04a")
            for mut, row in mutations:
                print(mut, row)

                row_path = pdb_path.pdb_thruputs + "back_" + str(row) + "/"
                print("### ... change directory", row_path)
                argus.params["thisrow"] = row
                argus.params["thismut"] = mut
                pdb_path.goto_job_dir(row_path, args, argus.params, "_inputs04a",emptyDirectory=True)
                print(
                    "### foldx03: ... copying file",
                    pdb_path.pdb_thruputs + pdbfile,
                    row_path + pdbfile,
                )
                copyfile(pdb_path.pdb_thruputs + pdbfile, row_path + pdbfile)

                fx_runner = Foldx.Foldx(argus.arg("foldxe"))                
                filename = pdb_path.pdb_inputs + "Coverage.csv"
                
                # prepare the coverage dataframe
                if exists(filename):
                    fdfp = FileDf.FileDf(filename)
                    cov_df = fdfp.openDataFrame()                    
                else:
                    empty_dic = {"source": [],"gene": [],"accession": [],"pdb": [],"method": [],"resolution": [],"chain": [],"pdb_start": [],"pdb_end": [],"gene_start": [],"gene_end": [],"coverage": [],"score": []}
                    cov_df = pd.DataFrame.from_dict(empty_dic)
                
                buildModel = True
                if buildModel: #run buildmodel

                    ddg_files = []
                    tag = 0
                    pdb_muts = mut.split(",")

                    allAtOnce=True
                    if allAtOnce:
                        muts = []
                        tag = 0
                        ddg_file = row_path + "Dif_" + str(tag) + "_" + pdbname + ".fxout"
                        for pm in pdb_muts:                            
                            muts.append(pm)                        
                        fx_runner.runBuild(pdbfile, pdb_muts, tag)                        
                        ddg_files.append([ddg_file, pdb_muts])
                    else:
                        for pm in pdb_muts:
                            tag += 1
                            ################################################################
                            fx_runner.runBuild(pdbfile, [pm], tag)
                            ################################################################
                            ddg_file = row_path + "Dif_" + str(tag) + "_" + pdbname + ".fxout"
                            ddg_files.append([ddg_file, pm])

                    # create them all into 1 ddg file
                    # df_file = row_path + "ddg_buildmodel.csv"
                    df_file = (
                        pdb_path.pdb_thruputs + "agg/" + str(row) + "_ddg_background.csv"
                    )
                    
                    fx_runner.createBuildCsv(
                        row_path, pdbname, pdb_muts, mut, cov_df, ddg_files, df_file, allAtOnce
                    )

                else:
                    ###########################################################################
                    fx_runner.runPosscan(pdbfile, mut)
                    ###########################################################################
                    pdb = pdbcode + "_rep" + str(argus.arg("repairs"))
                    # pass in the coverage to annotate the csv file
                    filename = pdb_path.pdb_inputs + "Coverage.csv"
                    ddg_file = row_path + "PS_" + pdb + "_scanning_output.txt"
                    df_file = (
                        pdb_path.pdb_thruputs + "agg/" + str(row) + "_ddg_background.csv"
                    )
                    fx_runner.createPosscanCsv(row_path, pdb, mut.split(","), [], cov_df, ddg_file, df_file)                    
        print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
