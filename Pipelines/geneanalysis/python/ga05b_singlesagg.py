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

# import from the shared library in Mutein/Pipelines/shared/lib
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/libs"
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

    work_path = pdb_path.pdb_thruputs + "vagg/"
    argus.params["work_path"] = work_path
    pdb_path.goto_job_dir(argus.arg("work_path"), args, argus.params, "_inputs07")
    ############################################
    params_file = pdb_path.pdb_thruputs + "singles_" + str(argus.arg("split")) + ".txt"
    fdfp = FileDf.FileDf(params_file, sep=" ", cols=["pdb", "gene_mut", "pdb_mut", "task"], header=False)
    pm_df = fdfp.openDataFrame()
    

                
    analyses = []
    analyses.append([   "ddg_posscan.csv",
                        pdb_path.pdb_outputs + "ddg_posscan.csv",                        
                        pdb_path.pdb_outputs + argus.arg("pdb")+"_ps_"+str(argus.arg("repairs"))+"_variants_plot.png"                        
                    ])
    analyses.append([   "ddg_buildmodel.csv",
                        pdb_path.pdb_outputs + "ddg_buildmodel.csv",                        
                        pdb_path.pdb_outputs + argus.arg("pdb")+"_bm_"+str(argus.arg("repairs"))+"_variants_plot.png"                        
                    ])
    
    for in_csv, out_csv, plot_file in analyses: 
        all_df = []
        for i in range(len(pm_df.index)):
            r = pm_df["task"][i]
            in_csv_i = work_path+str(r) + "_"+in_csv
            if exists(in_csv_i):
                fdf = FileDf.FileDf(in_csv_i)
                all_df.append(fdf.openDataFrame())
        if len(all_df) > 0:
            ddg_df = pd.concat(all_df, ignore_index=True)        
            ddg_df.to_csv(out_csv, index=False)        
            ana = Analysis.Analysis(ddg_df, argus.arg("pdb"))
            ana.createDdgResidue(plot_file, "variants", xax="gene_no")

        print("### COMPLETED FoldX aggregate job ###")
        print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline04"](sys.argv)
