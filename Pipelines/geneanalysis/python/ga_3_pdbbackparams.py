"""
-----------------------------
RSA 15/03/2022
-----------------------------
This file takes a pdb code (file must be located in the same directory in the format 1xyz.pdb)
It formats the pdb file into a parameter file suitable for foldx PositionScan
The output is in the same directory with the name
scanparams_1xyz.txt
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
# import sys
import os
import pandas as pd
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
def run_pipeline(args):
    print("### FoldX make params job ###")
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
        
    work_path = pdb_path.pdb_thruputs + "params" + str(argus.arg("split")) + "/"
    argus.addConfig({"work_path": work_path})
    pdb_path.goto_job_dir(argus.arg("work_path"), args, argus.params, "_inputs02")
    
    # chainid = argus.arg("chain")
    rows = int(argus.arg("split"))

    ##########################################################

    ##### Open the pdb file ################################
    pdb_file = pdbcode + ".pdb"        
    #pdb_file = pdbcode + "_rep" + str(argus.arg("repairs")) + ".pdb"        
    if exists(pdb_path.pdb_inputs+pdb_file):
        with open(pdb_path.pdb_inputs + pdb_file) as f:
            pdbcontent = f.readlines()

        ##### Amino acid dictionary to convert between 3 and 1 codes
        aa_dict = {
            "ALA": "A",
            "CYS": "C",
            "ASP": "D",
            "GLU": "E",
            "PHE": "F",
            "GLY": "G",
            "HIS": "H",
            "ILE": "I",
            "LYS": "K",
            "LEU": "L",
            "MET": "M",
            "ASN": "N",
            "PRO": "P",
            "GLN": "Q",
            "ARG": "R",
            "SER": "S",
            "THR": "T",
            "VAL": "V",
            "TRP": "W",
            "TYR": "Y",
        }

        ##### Go through each line and get the WT residue, then construct the foldx string to mutate it to everything, add to a list
        params_lst = []
        for line in pdbcontent:
            line = line.strip()
            if line.startswith("ATOM"):
                linecontents = [
                    line[:6],
                    line[6:11],
                    line[12:16],
                    line[17:20],
                    line[21],
                    line[22:26],
                    line[30:38],
                    line[38:46],
                    line[46:54],
                ]
                # this is e.g.   ATOM        4481         CB            THR         A         716
                #                          atom no        atom name     amino acid   chain   residue number
                atomtype = linecontents[2].strip()
                if atomtype == "CA":
                    chain = linecontents[4].strip()
                    if chain == chain:  # list of chains
                        aaa = linecontents[3].strip()
                        if aaa in aa_dict:
                            aa = aa_dict[aaa]
                            mut = (
                                aa + linecontents[4].strip() + linecontents[5].strip() + "a"
                            )  # aa chain rid mutation = mutation string
                            params_lst.append(mut)
                        else:
                            print("!Error maybe?", aaa)  # TODO think about this

        ##### Create a dataframe for the paramterfile in the number of chunks specified
        ##### Open up the coverage file
        #filename = pdb_path.pdb_inputs + "Coverage.csv"
        #fdfp = FileDf.FileDf(filename)
        #cov_df = fdfp.openDataFrame()
        #print(cov_df)

        rows = int(rows)
        param_dic = {}
        param_dic["pdb"] = []
        param_dic["mutation"] = []
        param_dic["row"] = []
        mut_str = ""

        total_muts = len(params_lst)
        chunk = int(total_muts / rows)
        remainder = int(total_muts % rows)
        # so until we get to the remainer we need chunk +1 on each row
        row_size = 0
        row = 0
        for i in range(len(params_lst)):
            mut = params_lst[i]
            if row_size == 0:
                param_dic["pdb"].append(pdbcode)
                param_dic["mutation"].append(mut)
                row += 1
                param_dic["row"].append("" + str(row))
            else:
                param_dic["mutation"][row - 1] = param_dic["mutation"][row - 1] + "," + mut
            row_size += 1

            if row_size == chunk and row > remainder:
                row_size = 0
            elif row_size == chunk + 1:
                row_size = 0
        print(total_muts, rows, chunk, row)
        ##### Turn the dictionary into a dataframe
        data_params = pd.DataFrame.from_dict(param_dic)
        filename = pdb_path.pdb_thruputs + "params_background.txt"
        print("### foldx02: ... savig df", filename)
        data_params.to_csv(filename, index=False, sep=" ", header=True)
        
    print("MUTEIN SCRIPT ENDED")
                    

##########################################################################################
if __name__ == "__main__":
    import sys
    globals()["run_pipeline"](sys.argv)
