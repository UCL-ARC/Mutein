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
import itertools

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
def run_pipeline05(args):
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
        work_path = pdb_path.pdb_thruputs + "vparams/"
        pdb_path.goto_job_dir(work_path, args, argus.params, "_inputs05")
        ############################################        
        variant = argus.arg("variant")
        # chainid = argus.arg("chain")
        splitrows = int(argus.arg("split"))

        # variant file is in the pdb inputs
        pdb_path = Paths.Paths(        
            data_dir,
            install_dir + "Pipelines/geneanalysis",
            dataset=dataset,
            gene=gene,
            pdb=pdbcode,
        )

        in_mutations_file = pdb_path.pdb_inputs + "variants.csv"
        new_mutations_file = pdb_path.pdb_outputs + "variants.csv"
        print("### foldx05: ... copying file", in_mutations_file, new_mutations_file)
        copyfile(in_mutations_file, new_mutations_file)
        ##### Open the variant file ################################
        variant_df = pd.read_csv(new_mutations_file)
        if variant != "*":
            mutations = variant_df.query("variant == '" + variant + "'")
        else:
            mutations = variant_df
        print(mutations)
        mut_list = []
        for i in range(len(mutations.index)):
            chain = mutations["chain"][i]
            mut = mutations["mutation"][i]
            pdb_mut = mutations["pdb_mut"][i]
            mut_list.append([mut, pdb_mut, chain])

        ##### Create a dataframe for the paramterfile in the number of chunks specified
        total_muts = len(mut_list)
        chunk = int(total_muts / splitrows)
        remainder = int(total_muts % splitrows)
        # so until we get to the remainer we need chunk +1 on each row
        param_dic = {}
        param_dic["pdb"] = []
        # param_dic["chain"] = []
        param_dic["gene_mut"] = []
        param_dic["pdb_mut"] = []
        param_dic["row"] = []
        row_size = 0
        row = 0
        for i in range(len(mut_list)):
            mut, pdb_mut, chain = mut_list[i]
            mutscan = mut[0] + chain + mut[1:]  # format for posscan
            pdb_mutscan = pdb_mut[0] + chain + pdb_mut[1:]  # format for posscan
            if row_size == 0:
                param_dic["pdb"].append(pdbcode)
                # param_dic["chain"].append(chainid)
                param_dic["gene_mut"].append(mutscan)
                param_dic["pdb_mut"].append(pdb_mutscan)
                row += 1
                param_dic["row"].append("" + str(row))
            else:
                param_dic["gene_mut"][row - 1] = (
                    param_dic["gene_mut"][row - 1] + "," + mutscan
                )
                param_dic["pdb_mut"][row - 1] = (
                    param_dic["pdb_mut"][row - 1] + "," + pdb_mutscan
                )
            row_size += 1

            if row_size == chunk and row > remainder:
                row_size = 0
            elif row_size == chunk + 1:
                row_size = 0
        print(total_muts, splitrows, chunk, row)
        ##### Turn the dictionary into a dataframe
        data_params = pd.DataFrame.from_dict(param_dic)
        filename = pdb_path.pdb_thruputs + "params_variants.txt"
        print("### foldx05: ... savig df", filename)
        data_params.to_csv(filename, index=False, sep=" ", header=True)
        all_params.append(data_params)
    
    all_path = gene_path.gene_outputs + "params_variants.txt"
    all_df = pd.concat(all_params, axis=0)
    all_df.to_csv(all_path,index=False,sep=" ",header=True)


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline05"](sys.argv)
