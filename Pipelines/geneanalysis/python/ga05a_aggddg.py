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
from os.path import exists
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

    gene_path = Paths.Paths(        
        data_dir,
        install_dir + "Pipelines/geneanalysis",
        dataset=dataset,
        gene=gene,        
    )
    pdbtasks = gene_path.gene_outputs + "pdb_tasklist.csv"
    fio = FileDf.FileDf(pdbtasks)
    df = fio.openDataFrame()
    
    for t in range(len(df.index)):
        pdbcode = df["pdb"][t].lower()
        
        pdb_path = Paths.Paths(        
            data_dir,
            install_dir + "Pipelines/geneanalysis",
            dataset=dataset,
            gene=gene,
            pdb=pdbcode,
        )
        # pdb_config = Config.Config(pdb_path.pdb_inputs + "/config.yml")
        # argus.addConfig(pdb_config.params)

        work_path = pdb_path.pdb_thruputs + "agg/"
        argus.params["work_path"] = work_path
        pdb_path.goto_job_dir(argus.arg("work_path"), args, argus.params, "_inputs05a")
        ############################################
        params_file = gene_path.gene_outputs +  "params_background.txt"
        fdfp = FileDf.FileDf(
            params_file, sep=" ", cols=["pdb", "mut", "task"], header=False
        )
        pm_df = fdfp.openDataFrame()
        all_df = []
        for i in range(len(pm_df.index)):
            r = pm_df["task"][i]
            rpdb = pm_df["pdb"][i]
            if rpdb == pdbcode:            
                # the file has already been turned into a dataframe called posscan_df.csv
                in_csv_i = work_path+str(r) + "_ddg_background.csv"                
                if exists(in_csv_i):
                    fdf = FileDf.FileDf(in_csv_i)
                    all_df.append(fdf.openDataFrame())
        if len(all_df) > 0:
            ddg_df = pd.concat(all_df, ignore_index=True)
            df_file = pdb_path.pdb_outputs + "ddg_background.csv"
            ddg_df.to_csv(df_file, index=False)

            plot_file = (
                pdb_path.pdb_outputs
                + argus.arg("pdb")
                + "_"
                + str(argus.arg("repairs"))
                + "_background_plot.png"
            )
            plot_file_gene = (
                pdb_path.pdb_outputs
                + argus.arg("pdb")
                + "_"
                + str(argus.arg("repairs"))
                + "_background_plot_gene.png"
            )

            ana = Analysis.Analysis(ddg_df, argus.arg("pdb"))
            ana.createDdgResidue(plot_file, "background")
            ana.createDdgResidue(
                plot_file_gene, "background muts", dropnagene=True, xax="gene_no"
            )

    print("### COMPLETED FoldX aggregate job ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline04"](sys.argv)
