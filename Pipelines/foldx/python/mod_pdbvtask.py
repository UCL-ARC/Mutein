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
import os
import pandas as pd
from shutil import copyfile
from os.path import exists
import sys

import _helper
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
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset","")
    gene = argus.arg("gene","")
    pdb = argus.arg("pdb", "")
    missing = argus.arg("missing", "N")
    repairs = str(argus.arg("repairs", "x"))

    task = int(argus.arg("task"))

    gene_path = Paths.Paths(
        data_dir,
        install_dir,
        dataset=dataset,
        gene=gene,
        pdb=pdb,
    )
    all_tasks = gene_path.inputs + "params_variants.txt"
    if missing.upper() == "Y":
        all_tasks = gene_path.inputs + "params_variants_incomplete.txt"
    fio = FileDf.FileDf(all_tasks, sep=" ", header=True)
    df = fio.openDataFrame()

    if task <= len(df.index):
        pdbcode = df["pdb"][task - 1].lower()

        pdb_path = Paths.Paths(
            data_dir,
            install_dir,
            dataset=dataset,
            gene=gene,
            pdb=pdbcode,
        )
        # pdb_config = Config.Config(pdb_path.pdb_inputs + "/config.yml")
        # argus.addConfig(pdb_config.params)
        task = argus.arg("task", "none")
        mutation_string = argus.arg("mutation", "none")
        ############################################
        pdbfile = pdbcode + "_rep" + repairs + ".pdb"
        mutations = []
        # task=all means all, task=1:n means an explicit row, row=-1 means the mutation string has been passd in explicitly
        if mutation_string == "none":
            if task == "all":
                for i in range(len(df.index)):
                    gene_mut = df["gene_mut"][i]
                    pdb_mut = df["pdb_mut"][i]
                    row = df["row"][i]
                    mutations.append([pdb_mut, gene_mut, row])
            else:
                if int(task) <= len(df.index):
                    gene_mut = df["gene_mut"][int(task) - 1]
                    pdb_mut = df["pdb_mut"][int(task) - 1]
                    row = df["row"][int(task) - 1]
                    mutations.append([pdb_mut, gene_mut, row])
        else:
            # we have specified a mutation and row from the file
            mutations.append([mutation_string, 0])

        work_path = pdb_path.pdb_thruputs + "vagg/"
        pdb_path.goto_job_dir(work_path, args, argus.params, "_inputs04b")

        for pdb_mut, gene_mut, row in mutations:
            print(pdb_mut, gene_mut, row)
            row_path = pdb_path.pdb_thruputs + "var_" + str(row) + "/"
            print("### ... change directory", row_path)
            argus.params["thisrow"] = row
            argus.params["genemut"] = gene_mut
            argus.params["pdbmut"] = pdb_mut
            pdb_path.goto_job_dir(row_path, args, argus.params, "_inputs06",emptyDirectory=True)
            print(
                "### foldx06: ... copying file",
                pdb_path.pdb_thruputs + pdbfile,
                row_path + pdbfile,
            )
            copyfile(pdb_path.pdb_inputs + pdbfile, row_path + pdbfile)

            #### TEMPORARILY DO BOTH POSCAN AND BUILD #####################
            fx_runner = Foldx.Foldx(argus.arg("foldxe"))
            pdb = pdbcode + "_rep" + repairs
            filename = pdb_path.pdb_inputs + "Coverage.csv"
            if exists(filename):
                fdfp = FileDf.FileDf(filename)
                cov_df = fdfp.openDataFrame()                    
            else:
                empty_dic = {"source": [],"gene": [],"accession": [],"pdb": [],"method": [],"resolution": [],"chain": [],"pdb_start": [],"pdb_end": [],"gene_start": [],"gene_end": [],"coverage": [],"score": []}
                cov_df = pd.DataFrame.from_dict(empty_dic)            
            if gene_mut.lower() == "x" or gene_mut == "":
                gene_muts = []
            else:
                gene_muts = gene_mut.split(",")
                        
            if False:
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~ POSITION SCAN ~~~~~~~~~~~~~~~~~~~~~~~~#
                ################################################################
                fx_runner.runPosscan(pdbfile, pdb_mut)
                ################################################################
                ddg_file = row_path + "PS_" + pdb + "_scanning_output.txt"
                # df_file = row_path + "ddg_posscan.csv"
                df_file = (
                    pdb_path.pdb_thruputs + "vagg/" + str(row) + "_ddg_posscan.csv"
                )
                pdb_muts = pdb_mut.split(",")
                fx_runner.createPosscanCsv(
                    row_path, pdb, pdb_muts, gene_muts, cov_df, ddg_file, df_file
                )

            if True:
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~ BUILD MODEL ~~~~~~~~~~~~~~~~~~~~~~~~#
                ddg_files = []
                tag = 0
                pdb_muts = pdb_mut.split(",")

                allAtOnce=True
                if allAtOnce:
                    muts = []
                    tag = 0
                    ddg_file = row_path + "Dif_" + str(tag) + "_" + pdb + ".fxout"
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
                        ddg_file = row_path + "Dif_" + str(tag) + "_" + pdb + ".fxout"
                        ddg_files.append([ddg_file, pm])


                # create them all into 1 ddg file
                # df_file = row_path + "ddg_buildmodel.csv"
                df_file = (
                    pdb_path.pdb_thruputs + "vagg/" + str(row) + "_ddg_buildmodel.csv"
                )
                fx_runner.createBuildCsv(
                    row_path, pdb, pdb_muts, gene_muts, cov_df, ddg_files, df_file,allAtOnce
                )
            import mod_cleaning as cleaner        
            cleaner.deletePath(row_path,remPath=True)
    
    #print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
