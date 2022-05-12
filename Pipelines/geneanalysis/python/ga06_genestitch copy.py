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
    gene_path = Paths.Paths(        
        data_dir,
        install_dir + "Pipelines/geneanalysis",
        dataset=dataset,
        gene=gene,
    )
    work_path = gene_path.gene_outputs + "agg/"
    argus.params["work_path"] = work_path
    gene_path.goto_job_dir(argus.arg("work_path"), args, argus.params, "_inputs01")
    ############################################
    params_file = gene_path.gene_outputs + "Coverage_all.csv"
    fdfp = FileDf.FileDf(params_file)    
    pm_df = fdfp.openDataFrame()    
    pdbs = pm_df["pdb"].unique()
    all_df_vars_ps = []
    all_df_vars_bm = []
    all_df_backs_all = []
    all_df_backs = []
    all_df_vars_AF_ps = []
    all_df_vars_AF_bm = []
    all_df_backs_AF = []
    for pdb in pdbs:
        r = pdb.lower()
        # the file has already been turned into a dataframe called posscan_df.csv
        pdb_path = Paths.Paths(            
            data_dir,
            install_dir + "Pipelines/geneanalysis",
            dataset=dataset,
            gene=gene,
            pdb=r,
        )
        file_var_ps = pdb_path.pdb_outputs + "ddg_posscan.csv"
        file_var_bm = pdb_path.pdb_outputs + "ddg_buildmodel.csv"
        file_back = pdb_path.pdb_outputs + "ddg_background.csv"
        if exists(file_var_ps):
            fdf = FileDf.FileDf(file_var_ps)            
            if r[:2].upper() == "AF":
                all_df_vars_AF_ps.append(fdf.openDataFrame())
            else:
                all_df_vars_ps.append(fdf.openDataFrame())
        if exists(file_var_bm):
            fdf = FileDf.FileDf(file_var_bm)            
            if r[:2].upper() == "AF":
                all_df_vars_AF_bm.append(fdf.openDataFrame())
            else:
                all_df_vars_bm.append(fdf.openDataFrame())
        if exists(file_back):
            fdf = FileDf.FileDf(file_back)
            if r[:2].upper() == "AF":
                all_df_backs_AF.append(fdf.openDataFrame())
            else:
                all_df_backs.append(fdf.openDataFrame())
            all_df_backs_all.append(fdf.openDataFrame())

    # NON Alpha-Fold structures
    outputs = []
    outputs.append(#variants done with posscan
        [
            all_df_vars_ps,
            "ddg_variants_ps.csv",
            "ddg_variants_ps.png",
            "variants",
            "gene_no",
            False,            
        ]
    )
    outputs.append(#variants done with buildmodel
        [
            all_df_vars_bm,
            "ddg_variants_bm.csv",
            "ddg_variants_bm.png",
            "variants",
            "gene_no",
            False,            
        ]
    )
    outputs.append(
        [
            all_df_backs,
            "ddg_background.csv",
            "ddg_background.png",
            "background",
            "pdb_rid",
            False,            
        ]
    )
    outputs.append(
        [
            all_df_backs,
            "ddg_background.csv",
            "ddg_background_gene.png",
            "background gene",
            "gene_no",
            True,            
        ]
    )
    outputs.append(
        [
            all_df_vars_AF_ps,
            "AF_ddg_variants_ps.csv",
            "AF_ddg_variants_ps.png",
            "variants AF",
            "gene_no",
            False,            
        ]
    )
    outputs.append(
        [
            all_df_vars_AF_bm,
            "AF_ddg_variants_bm.csv",
            "AF_ddg_variants_bm.png",
            "variants AF",
            "gene_no",
            False,            
        ]
    )
    outputs.append(
        [
            all_df_backs_AF,
            "AF_ddg_background.csv",
            "AF_ddg_background.png",
            "background AF",
            "pdb_rid",
            False,            
        ]
    )
    outputs.append(
        [
            all_df_backs_AF,
            "AF_ddg_background.csv",
            "AF_ddg_background_gene.png",
            "background AF gene",
            "gene_no",
            True,            
        ]
    )

    for dfs, csv_file, png_file, title, xax, nagene in outputs:
        if len(dfs) > 0:
            ddg_df = pd.concat(dfs, ignore_index=True)
            df_file = gene_path.gene_outputs + csv_file
            ddg_df.to_csv(df_file, index=False)
            plot_file = gene_path.gene_outputs + png_file
            anav = Analysis.Analysis(ddg_df, gene)
            anav.createDdgResidue(plot_file, title, xax=xax, dropnagene=nagene)
    # The coverage of the pdb structures is across all structures, inc AF and SM
    if len(all_df_backs_all) > 0:
            ddg_df = pd.concat(all_df_backs_all, ignore_index=True)
            df_file = gene_path.gene_outputs + "all_background.csv"
            ddg_df.to_csv(df_file, index=False)
            plot_file = gene_path.gene_outputs + "all_pdbcoverage.png"
            anav = Analysis.Analysis(ddg_df, gene)                        
            anav.createPdbSummary(plot_file, "PDB Coverage")

    print("### COMPLETED GENE STITCH job ###")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
