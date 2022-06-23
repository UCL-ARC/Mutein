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


def run_pipeline(args):
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
    pdbcode = argus.arg("pdb", "").lower()

    gene_path = Paths.Paths(
        data_dir,
        install_dir + "Pipelines/geneanalysis",
        dataset=dataset,
        gene=gene,
        pdb=pdbcode,
    )

    pdb_list = []
    if pdbcode != "":
        pdb_list.append(pdbcode)
    else:
        pdbtasks = gene_path.gene_outputs + "pdb_tasklist.csv"
        fio = FileDf.FileDf(pdbtasks)
        df = fio.openDataFrame()

        for t in range(len(df.index)):
            pdbcode = df["pdb"][t].lower()
            pdb_list.append(pdbcode)

    for pdbcode in pdb_list:

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
        pdb_path.goto_job_dir(argus.arg("work_path"), args, argus.params, "_inputs05b")
        ############################################
        params_file = gene_path.thruputs + "params_variants.txt"
        fdfp = FileDf.FileDf(params_file, sep=" ", header=True)
        pm_df = fdfp.openDataFrame()

        analyses = []
        analyses.append(
            [
                "ddg_posscan.csv",
                pdb_path.pdb_outputs + "ddg_posscan.csv",
                pdb_path.pdb_outputs
                + argus.arg("pdb")
                + "_ps_"
                + str(argus.arg("repairs"))
                + "_variants_plot.png",
            ]
        )
        analyses.append(
            [
                "ddg_buildmodel.csv",
                pdb_path.pdb_outputs + "ddg_buildmodel.csv",
                pdb_path.pdb_outputs
                + argus.arg("pdb")
                + "_bm_"
                + str(argus.arg("repairs"))
                + "_variants_plot.png",
            ]
        )

        for in_csv, out_csv, plot_file in analyses:
            all_df = []
            exists_all = True
            for i in range(len(pm_df.index)):
                r = pm_df["row"][i]
                rpdb = pm_df["pdb"][i]
                if rpdb == pdbcode:
                    in_csv_i = work_path + str(r) + "_" + in_csv
                    if exists(in_csv_i):
                        fdf = FileDf.FileDf(in_csv_i)
                        all_df.append(fdf.openDataFrame())
                    else:
                        exists_all = False
            if len(all_df) > 0 and exists_all:
                ddg_df = pd.concat(all_df, ignore_index=True)
                ddg_df.to_csv(out_csv, index=False)
                ana = Analysis.Analysis(ddg_df, argus.arg("pdb"))
                ana.createDdgResidue(plot_file, "variants", xax="gene_no")
            else:
                print("vAgg - not all results present for aggregation",in_csv)

    print("### COMPLETED FoldX aggregate job ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
