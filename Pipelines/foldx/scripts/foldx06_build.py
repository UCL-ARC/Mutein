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


#import from the shared library in Mutein/Pipelines/shared/lib
import sys
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + '/shared/libs'
sys.path.append(retpath)
import Paths
import Arguments
import Config
import Foldx


def run_pipeline06(args):
    print("### FoldX build job ###")
    print(args)    
    argus = Arguments.Arguments(args)
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    pdbcode = argus.arg("pdb").lower()
    pdb_path = Paths.Paths("pdb",dataset=dataset,gene=gene,pdb=pdbcode)
    #pdb_config = Config.Config(pdb_path.pdb_inputs + "/config.yml")
    #argus.addConfig(pdb_config.params)    
    task = argus.arg("task","none")
    mutation_string = argus.arg("mutation","none")

    ############################################
    # set up the files and directories
    pdbfile = pdbcode + "_rep" + str(argus.arg("repairs")) + ".pdb"
    mutations = []

    if mutation_string == "none":
        filename = pdb_path.pdb_thruputs + "singles_" + str(argus.arg("split")) + ".txt"
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

        row_path = pdb_path.pdb_thruputs + str(argus.arg("split")) + "_" + str(row) + "_build" + "/"
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
        fx_runner.runBuild(pdbfile,mut,15)    
        ddg_file = row_path+"/Dif_" + pdbfile + ".fxout"
        df_file = row_path+"/pdbfile_build_DDG.csv"        
        fx_runner.createBuildCsv(pdbcode,mut,ddg_file,df_file,task)

if __name__ == "__main__":
    import sys

    globals()["run_pipeline06"](sys.argv)
