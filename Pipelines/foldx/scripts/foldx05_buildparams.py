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

#import from the shared library in Mutein/Pipelines/shared/lib
import sys
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + '/shared/libs'
sys.path.append(retpath)
import Paths
import Arguments
import Config


##### INPUTS #############################################
# The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
def run_pipeline05(args):
    print("### FoldX make variant params job ###")
    print(args)
    ##############################################    
    argus = Arguments.Arguments(args)
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    pdbcode = argus.arg("pdb")
    pdb_path = Paths.Paths("pdb",dataset=dataset,gene=gene,pdb=pdbcode)
    work_path = pdb_path.pdb_thruputs + "vparams/"
    pdb_path.goto_job_dir(work_path, args, argus.params, "_inputs05")    
    ############################################
    pdb = argus.arg("pdb")
    #jobname = argus.arg("name")
    #row = argus.arg("row","row0")
    variant = argus.arg("variant")
    chainid = argus.arg("chain")
    
    # variant file is in the pdb inputs
    pdb_path = Paths.Paths("pdb",dataset=dataset,gene=gene,pdb=pdbcode)

    in_mutations_file = pdb_path.pdb_inputs + "variants.csv"
    new_mutations_file = pdb_path.pdb_outputs + "variants.csv"
    print("### foldx05: ... copying file", in_mutations_file, new_mutations_file)
    copyfile(in_mutations_file, new_mutations_file)
    ##### Open the variant file ################################
    variant_df = pd.read_csv(new_mutations_file)
    if variant !="*":
        mutations = variant_df.query("variant == '" + variant + "'")
    else:
        mutations = variant_df
    print(mutations)
    mut_list = []
    for i in range(len(variant_df.index)):
        mut = variant_df["mutation"][i]
        mut_list.append(mut)

    mut_combos = []
    for i in range(len(mut_list)):
        for subset in itertools.combinations(mut_list, i + 1):
            mut_combos.append(subset)
            print(subset)

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
        print(pdb,chainid,mut_str,"row" + str(len(mut)) + "_var" + str(row))

    ##### Turn the dictionary into a dataframe
    data_params = pd.DataFrame.from_dict(param_dic)
    filename = pdb_path.pdb_thruputs + "variant_params.txt"
    print("## Foldx.vparam: exporting variant file to",filename)
    data_params.to_csv(filename, index=False, sep=" ", header=False)


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline05"](sys.argv)
