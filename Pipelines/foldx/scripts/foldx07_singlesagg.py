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
from os.path import exists

#import from the shared library in Mutein/Pipelines/shared/lib
import sys
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + '/shared/libs'
sys.path.append(retpath)
import Paths
import Arguments
import Config
import Analysis
import FileDf



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

    fdfp = FileDf.FileDf(params_file,sep=" ",cols=["pdb","gene_mut","pdb_mut","task"],header=False)
    pm_df = fdfp.openDataFrame()
    all_df = []
    for i in range(len(pm_df.index)):
        r = pm_df["task"][i]
        # the file has already been turned into a dataframe called posscan_df.csv        
        file_path = pdb_path.pdb_thruputs + str(argus.arg("split")) + "_" + str(r) + "_var/posscan_df.csv" 
        if exists(file_path):
            fdf = FileDf.FileDf(file_path)
            all_df.append(fdf.openDataFrame())        
    ddg_df = pd.concat(all_df,ignore_index=True)        
    df_file = (
        pdb_path.pdb_outputs + "ddg_variants.csv"
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
    ana.createDdgResidue(plot_file,"variants")
    
    print("### COMPLETED FoldX aggregate job ###")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline04"](sys.argv)
