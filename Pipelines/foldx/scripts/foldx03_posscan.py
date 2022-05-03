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

#import from the shared library in Mutein/Pipelines/shared/lib
import sys
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + '/shared/libs'
sys.path.append(retpath)
import Paths
import Arguments
import Config
import Foldx


##### INPUTS #############################################
# The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)


def run_pipeline03(args):
    print("### FoldX position scan job ###")
    print(args)
    argus = Arguments.Arguments(args)
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    pdbcode = argus.arg("pdb")
    pdb_path = Paths.Paths("pdb",dataset=dataset,gene=gene,pdb=pdbcode)
    #pdb_config = Config.Config(pdb_path.pdb_inputs + "/config.yml")
    #argus.addConfig(pdb_config.params)    
    task = argus.arg("task","none")
    mutation_string = argus.arg("mutation","none")
    ############################################
    pdbfile = pdbcode + "_rep" + str(argus.arg("repairs")) + ".pdb"
    mutations = []
    # task=all means all, task=1:n means an explicit row, row=-1 means the mutation string has been passd in explicitly
    if mutation_string == "none":
        filename = pdb_path.pdb_thruputs + "params_" + str(argus.arg("split")) + ".txt"
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
                if int(task) <= len(paramscontent):
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

        row_path = pdb_path.pdb_thruputs + str(argus.arg("split")) + "_" + str(row) + "/"
        print("### ... change directory", row_path)
        argus.params["thisrow"] = row
        argus.params["thismut"] = mut
        pdb_path.goto_job_dir(row_path, args, argus.params, "_inputs03")
        print(
            "### foldx03: ... copying file",
            pdb_path.pdb_thruputs + pdbfile,
            row_path + pdbfile,
        )
        copyfile(pdb_path.pdb_thruputs + pdbfile, row_path + pdbfile)

        fx_runner = Foldx.Foldx(argus.arg("foldxe"))    
        fx_runner.runPosscan(pdbfile,mut)

        


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline03"](sys.argv)
