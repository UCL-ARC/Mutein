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
    
    all_back = []
    all_var_scan = []
    all_var_build = []
    
    # first concat all the dfs, variants and backgrount
    for pdb in pdbs:        
        # the file has already been turned into a dataframe called posscan_df.csv
        pdb_path = Paths.Paths(            
            data_dir,
            install_dir + "Pipelines/geneanalysis",
            dataset=dataset,
            gene=gene,
            pdb=pdb,
        )
        file_var_ps = pdb_path.pdb_outputs + "ddg_posscan.csv"
        file_var_bm = pdb_path.pdb_outputs + "ddg_buildmodel.csv"
        file_back = pdb_path.pdb_outputs + "ddg_background.csv"
        if exists(file_var_ps):
            fdf = FileDf.FileDf(file_var_ps)                        
            all_var_scan.append(fdf.openDataFrame())            
        if exists(file_var_bm):
            fdf = FileDf.FileDf(file_var_bm)                        
            all_var_build.append(fdf.openDataFrame())            
        if exists(file_back):
            fdf = FileDf.FileDf(file_back)            
            all_back.append(fdf.openDataFrame())

    ddg_df_back = pd.concat(all_back,ignore_index=True)
    ddg_df_var_scan = pd.concat(all_var_scan,ignore_index=True)
    ddg_df_var_build = pd.concat(all_var_build,ignore_index=True)
    ddg_df_back.to_csv(gene_path.gene_outputs+"ddg_background.csv", index=False)
    ddg_df_var_scan.to_csv(gene_path.gene_outputs+"ddg_variant_ps.csv", index=False)
    ddg_df_var_build.to_csv(gene_path.gene_outputs+"ddg_variant_bm.csv", index=False)
            
    #DB Coverage reports    
    plot_file = gene_path.gene_outputs + "ALL_coverage.png"
    anav = Analysis.Analysis(ddg_df_back, gene)                        
    anav.createPdbSummary(plot_file, "PDB Coverage")
        
    tags = ["SMHOM","SMEXP","AF","EXP"]
    for tag in tags:
        df_back = ddg_df_back.query("source=='"+tag+"'")
        df_scan = ddg_df_var_scan.query("source=='"+tag+"'")
        df_build = ddg_df_var_build.query("source=='"+tag+"'")        
        
        anav_back = Analysis.Analysis(df_back, gene)
        anav_back.createDdgResidue(gene_path.gene_outputs+tag+ "_background_rid.png", "Background (residue) "+ tag, xax="pdb_rid", dropnagene=False)
        anav_back.createDdgResidue(gene_path.gene_outputs+tag+ "_background_gene.png", "Background (gene no) "+ tag, xax="gene_no", dropnagene=True)

        anav_scan = Analysis.Analysis(df_scan, gene)
        anav_scan.createDdgResidue(gene_path.gene_outputs+tag+ "_variant_ps.png", "Variant PosScan "+ tag, xax="gene_no", dropnagene=False)

        anav_build = Analysis.Analysis(df_build, gene)
        anav_build.createDdgResidue(gene_path.gene_outputs+tag+ "_variant_bm.png", "Variant BuildModel "+ tag, xax="gene_no", dropnagene=False)
                            
    print("### COMPLETED GENE STITCH job ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
