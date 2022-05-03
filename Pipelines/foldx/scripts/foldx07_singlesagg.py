"""
-----------------------------
RSA 17/03/22
-----------------------------

This aggregates the outputs from split positionscans into 1 file
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
import Analysis



def run_pipeline04(args):
    print("### Foldx aggregate pos scan ###")
    print(args)
    ##############################################
    argus = Arguments.Arguments(args)
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    pdbcode = argus.arg("pdb")
    pdb_path = Paths.Paths("pdb",dataset=dataset,gene=gene,pdb=pdbcode)
    #pdb_config = Config.Config(pdb_path.pdb_inputs + "/config.yml")
    #argus.addConfig(pdb_config.params) 

    work_path = pdb_path.pdb_thruputs + "agg/"
    argus.params["work_path"] = work_path
    pdb_path.goto_job_dir(argus.arg("work_path"), args, argus.params, "_inputs07")
    ############################################
    params_file = pdb_path.pdb_thruputs + "singles_" + str(argus.arg("split"))+ ".txt"
    rownum = 1
    with open(params_file) as fr:
        paramscontent = fr.readlines()
        rownum = len(paramscontent)

    ddg_file = (
        "PS_"
        + argus.arg("pdb")
        + "_rep"
        + str(argus.arg("repairs"))
        + "_scanning_output.txt"
    )
    with open(ddg_file, "w") as fw:
        for r in range(rownum):
            jobresults_file = (
                pdb_path.pdb_thruputs + str(argus.arg("split")) + "_" + str(r + 1) + "_var/" + ddg_file
            )
            if os.path.exists(jobresults_file):
                print("### foldx.aggregate",jobresults_file)
                with open(jobresults_file) as fr:
                    jobcontent = fr.readlines()
                fw.writelines(jobcontent)
            #else:
            #    print("No file",jobresults_file)

    # Make a dataframe
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
        "H1S": "o",
        "H2S": "e",
    }
    a_dict = {
        "A": "ALA",
        "C": "CYS",
        "D": "ASP",
        "E": "GLU",
        "F": "PHE",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "K": "LYS",
        "L": "LEU",
        "M": "MET",
        "N": "ASN",
        "P": "PRO",
        "Q": "GLN",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "V": "VAL",
        "W": "TRP",
        "Y": "TYR",
        "o": "H1S",
        "e": "H2S",
    }
    ddg_dic = {"aa": [], "chain": [], "rid": [], "mut": [], "mutid": [], "ddg": []}
    with open(ddg_file, "r") as fr:
        lines = fr.readlines()
    for line in lines:
        if len(line) > 7:
            lns = line.strip().split("\t")
            res = lns[0]
            aa = res[:3]
            chain = res[3:4]
            rid = res[4:-1]
            mut = res[-1:]
            ddg = lns[1]
            # print(aa,chain,rid,mut,ddg)
            ddg_dic["aa"].append(aa)
            ddg_dic["chain"].append(chain)
            ddg_dic["rid"].append(int(rid))
            ddg_dic["mut"].append(mut)
            ddg_dic["mutid"].append(aa_dict[aa] + chain + rid + mut)
            ddg_dic["ddg"].append(float(ddg))
    import pandas as pd

    ddg_df = pd.DataFrame.from_dict(ddg_dic)
    df_file = (
        pdb_path.pdb_outputs
        + argus.arg("pdb")
        + "_"
        + str(argus.arg("repairs"))
        + "_ddg_dataframe_var.csv"
    )
    ddg_df.to_csv(df_file, index=False)        
    plot_file = (
        pdb_path.pdb_outputs
        + argus.arg("pdb")
        + "_"
        + str(argus.arg("repairs"))
        + "_variants_plot.png"
    )
    ana = Analysis.Analysis(ddg_df,argus.arg("pdb"))
    ana.createDdgResidue(plot_file)
    
    print("### COMPLETED FoldX aggregate job ###")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline04"](sys.argv)
