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
import helper as hlp
import Arguments

##### INPUTS #############################################
# The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
def run_pipeline05(args):
    print("### FoldX make variant params job ###")
    print(args)
    ##############################################
    argus = Arguments.Arguments(args)
    work_path = argus.params["interim_path"] + "vparams/"
    argus.params["work_path"] = work_path
    hlp.goto_job_dir(argus.arg("work_path"), args, argus.params, "_inputs05")
    ############################################
    pdb = argus.arg("pdb")
    jobname = argus.arg("name")
    row = argus.arg("row")
    variant = argus.arg("variant")
    chainid = argus.arg("chain")
    # we want to work in the node directory first, the main pdb input file is a 1-off and lives in github (at the moment)
    in_mutations_file = argus.arg("input_path") + argus.arg("variantfile") + ".txt"
    new_mutations_file = argus.arg("output_path") + argus.arg("variantfile") + ".txt"
    print("### ... copying file", in_mutations_file, new_mutations_file)
    copyfile(in_mutations_file, new_mutations_file)
    ##### Open the variant file ################################
    variant_df = pd.read_csv(new_mutations_file)
    mutations = variant_df.query("Variant == '" + variant + "'")
    print(mutations)
    mut_list = []
    for i in range(len(variant_df.index)):
        mut = variant_df["Mutation"][i]
        mut_list.append(mut)

    mut_combos = []
    for i in range(len(mut_list)):
        for subset in itertools.combinations(mut_list, i + 1):
            mut_combos.append(subset)

    ##### Create a dataframe for the paramterfile in the number of chunks specified
    param_dic = {}
    param_dic["pdb"] = []
    param_dic["chain"] = []
    param_dic["mutation"] = []
    param_dic["row"] = []
    mut_str = ""
    count = -1
    row = 0
    for mut in mut_combos:
        row += 1
        mut_str = ""
        for m in mut:
            if len(mut_str) > 0:
                mut_str += ","
            mut_str += m
        param_dic["pdb"].append(pdb)
        param_dic["chain"].append(chainid)
        param_dic["mutation"].append(mut_str + ";")
        param_dic["row"].append("row" + str(len(mut)) + "_var" + str(row))

    ##### Turn the dictionary into a dataframe
    data_params = pd.DataFrame.from_dict(param_dic)
    filename = argus.arg("thruput_path") + "variant_params.txt"
    data_params.to_csv(filename, index=False, sep=" ", header=False)


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline05"](sys.argv)
