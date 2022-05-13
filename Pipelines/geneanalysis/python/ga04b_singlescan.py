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


def run_pipeline03(args):
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
    pdbcode = argus.arg("pdb").lower()
    pdb_path = Paths.Paths(        
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
    pdbfile = pdbcode + "_rep" + str(argus.arg("repairs")) + ".pdb"
    mutations = []
    # task=all means all, task=1:n means an explicit row, row=-1 means the mutation string has been passd in explicitly
    if mutation_string == "none":
        filename = pdb_path.pdb_thruputs + "singles_" + str(argus.arg("split")) + ".txt"
        fio = FileDf.FileDf(
            filename, sep=" ", cols=["pdb", "gene_mut", "pdb_mut", "task"], header=False
        )
        df = fio.openDataFrame()
        if task == "all":
            for i in range(len(df.index)):
                gene_mut = df["gene_mut"][i]
                pdb_mut = df["pdb_mut"][i]
                row = df["task"][i]
                mutations.append([pdb_mut, gene_mut, row])
        else:
            if int(task) <= len(df.index):
                gene_mut = df["gene_mut"][int(task) - 1]
                pdb_mut = df["pdb_mut"][int(task) - 1]
                row = df["task"][int(task) - 1]
                mutations.append([pdb_mut, gene_mut, row])
    else:
        # we have specified a mutation and row from the file
        mutations.append([mutation_string, 0])

    for pdb_mut, gene_mut, row in mutations:
        print(pdb_mut, gene_mut, row)
        row_path = (
            pdb_path.pdb_thruputs
            + str(argus.arg("split"))
            + "_"
            + str(row)
            + "_var"
            + "/"
        )
        print("### ... change directory", row_path)
        argus.params["thisrow"] = row
        argus.params["genemut"] = gene_mut
        argus.params["pdbmut"] = pdb_mut
        pdb_path.goto_job_dir(row_path, args, argus.params, "_inputs06")
        print(
            "### foldx06: ... copying file",
            pdb_path.pdb_thruputs + pdbfile,
            row_path + pdbfile,
        )
        copyfile(pdb_path.pdb_thruputs + pdbfile, row_path + pdbfile)

        #### TEMPORARILY DO BOTH POSCAN AND BUILD #####################
        fx_runner = Foldx.Foldx(argus.arg("foldxe"))
        pdb = pdbcode + "_rep" + str(argus.arg("repairs"))
        filename = pdb_path.pdb_inputs + "Coverage.csv"
        fdfp = FileDf.FileDf(filename)
        cov_df = fdfp.openDataFrame()        
        if gene_mut.lower() == "x" or gene_mut == "":
            gene_muts = []
        else:
            gene_muts = gene_mut.split(",")
        
        work_path = pdb_path.pdb_thruputs + "vagg/"
        pdb_path.goto_job_dir(work_path, args, argus.params, "_inputs04b")
        
        if True:
            #~~~~~~~~~~~~~~~~~~~~~~~~~~ POSITION SCAN ~~~~~~~~~~~~~~~~~~~~~~~~#        
            ################################################################
            #fx_runner.runPosscan(pdbfile, pdb_mut)
            ################################################################                        
            ddg_file = row_path + "PS_" + pdb + "_scanning_output.txt"
            #df_file = row_path + "ddg_posscan.csv"        
            df_file = pdb_path.pdb_thruputs+"vagg/"+str(row) + "_ddg_posscan.csv"              
            pdb_muts = pdb_mut.split(",")
            fx_runner.createPosscanCsv(row_path, pdb, pdb_muts, gene_muts, cov_df, ddg_file,df_file)
        
        if True:            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~ BUILD MODEL ~~~~~~~~~~~~~~~~~~~~~~~~#                                            
            ddg_files = []
            tag = 0
            pdb_muts = pdb_mut.split(",")
            for pm in pdb_muts:
                tag+=1
                ################################################################        
                #fx_runner.runBuild(pdbfile, pm,tag)
                ################################################################                        
                ddg_file = row_path + "Dif_" + str(tag) + "_" + pdb + ".fxout"                
                ddg_files.append([ddg_file,pm])
            
            #create them all into 1 ddg file
            #df_file = row_path + "ddg_buildmodel.csv"  
            df_file = pdb_path.pdb_thruputs+"vagg/"+str(row) + "_ddg_buildmodel.csv"              
            fx_runner.createBuildCsv(row_path, pdb, pdb_muts, gene_muts, cov_df, ddg_files,df_file)


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline03"](sys.argv)
