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
import sys

import _helper
import Paths
import Arguments
import Config
import Analysis
import FileDf
import ScoringMetric


def run_pipeline(args):
    print("### Foldx aggregate pos scan ###")
    print(args)
    ##############################################
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    gene_path = Paths.Paths(
        data_dir,
        install_dir,
        dataset=dataset,
        gene=gene,
    )
    work_path = gene_path.gene_outputs + "agg/"
    argus.params["work_path"] = work_path
    gene_path.goto_job_dir(argus.arg("work_path"), args, argus.params, "_inputs01")
    ############################################
    #params_file = gene_path.gene_outputs + "Coverage_all.csv"
    #fdfp = FileDf.FileDf(params_file)
    #pm_df = fdfp.openDataFrame()
    
    pdbs = []
    pdbtasks = gene_path.gene_outputs + "pdb_tasklist.csv"
    fio = FileDf.FileDf(pdbtasks)
    df = fio.openDataFrame()
    for t in range(len(df.index)):
        pdbcode = df["pdb"][t].lower()
        if pdbcode[0] != "#":
            pdbs.append(pdbcode)
    
    all_back_bm = []    
    all_var_build = []

    # first concat all the dfs, variants and backgrount
    exists_all_back = True
    exists_all_var = True
    for pdb in pdbs:
        # the file has already been turned into a dataframe called posscan_df.csv
        pdb_path = Paths.Paths(
            data_dir,
            install_dir,
            dataset=dataset,
            gene=gene,
            pdb=pdb,
        )        
        
        
        file_lst_var = pdb_path.pdb_inputs + "params_variants.txt"
        file_var_bm = pdb_path.pdb_outputs + "ddg_variants.csv"        
        file_back_bm = pdb_path.pdb_outputs + "ddg_background.csv"
        
        if exists(file_lst_var):
            fdf = FileDf.FileDf(file_lst_var)
            lstvardf = fdf.openDataFrame()
            if len(lstvardf.index) > 0:                 
                if exists(file_var_bm):
                    fdf = FileDf.FileDf(file_var_bm)                    
                    all_var_build.append(fdf.openDataFrame())                                    
                else:
                    exists_all_var = False                    
            else:                
                print("No variants to aggregate",gene,pdb)                    
        else:
            exists_all_var = False
        if exists(file_back_bm):
            fdf = FileDf.FileDf(file_back_bm)
            all_back_bm.append(fdf.openDataFrame())
        else:
            print("No pdb file for",file_back_bm)
            exists_all_back = False
        

    metric = ScoringMetric.ScoringMetric(gene_path,dataset,gene)
    if exists_all_back and len(all_back_bm) > 0:        
        ddg_df_back_bm = pd.concat(all_back_bm, ignore_index=True)
        ddg_df_back_bm['pdb'] = ddg_df_back_bm.apply(lambda row: metric.cutPdb(row['pdb']), axis=1)
        ddg_df_back_bm['score'] = ddg_df_back_bm.apply(lambda row: metric.getScore(row['pdb'])[0], axis=1)
        ddg_df_back_bm['method'] = ddg_df_back_bm.apply(lambda row: metric.getScore(row['pdb'])[1], axis=1)
        ddg_df_back_bm['resolution'] = ddg_df_back_bm.apply(lambda row: metric.getScore(row['pdb'])[2], axis=1)
        ddg_df_back_bm['coverage'] = ddg_df_back_bm.apply(lambda row: metric.getScore(row['pdb'])[3], axis=1)
        ddg_df_back_bm.to_csv(gene_path.gene_outputs + "ddg_background.csv", index=False)
        # DB Coverage reports
        #plot_file = gene_path.gene_outputs + "ALL_coverage.png"
        #anav = Analysis.Analysis(ddg_df_back_bm, gene)
        #anav.createPdbSummary(plot_file, "PDB Coverage")
    else:
        print("Gene stitch: not all pdb background files present for aggregation")

    if exists_all_var:
        if len(all_var_build) > 0:        
            ddg_df_var_build = pd.concat(all_var_build, ignore_index=True)
            ddg_df_var_build['pdb'] = ddg_df_var_build.apply(lambda row: metric.cutPdb(row['pdb']), axis=1)
            ddg_df_var_build['score'] = ddg_df_var_build.apply(lambda row: metric.getScore(row['pdb'])[0], axis=1)
            ddg_df_var_build['method'] = ddg_df_var_build.apply(lambda row: metric.getScore(row['pdb'])[1], axis=1)
            ddg_df_var_build['resolution'] = ddg_df_var_build.apply(lambda row: metric.getScore(row['pdb'])[2], axis=1)
            ddg_df_var_build['coverage'] = ddg_df_var_build.apply(lambda row: metric.getScore(row['pdb'])[3], axis=1)
            ddg_df_var_build.to_csv(
                gene_path.gene_outputs + "ddg_variants.csv", index=False
            )
        else:
            print("Gene stitch: no files available for aggreation")
    else:
        print("Gene stitch: not all pdb variant files present for aggregation")

    """
    tags = ["SMHOM", "SMEXP", "AF", "EXP"]
    for tag in tags:                        
        if exists_all_back:
            df_back = ddg_df_back_bm.query("source=='" + tag + "'")        
            anav_back = Analysis.Analysis(df_back, gene)
            anav_back.createDdgResidue(
                gene_path.gene_outputs + tag + "_background_rid.png",
                "Background (residue) " + tag,
                xax="pdb_rid",
                dropnagene=False,
            )
            anav_back.createDdgResidue(
                gene_path.gene_outputs + tag + "_background_gene.png",
                "Background (gene no) " + tag,
                xax="gene_no",
                dropnagene=True,
            )

        if exists_all_var and len(all_var_build) > 0:
            df_build = ddg_df_var_build.query("source=='" + tag + "'")            
            
            anav_build = Analysis.Analysis(df_build, gene)
            anav_build.createDdgResidue(
                gene_path.gene_outputs + tag + "_variant_bm.png",
                "Variant BuildModel " + tag,
                xax="gene_no",
                dropnagene=False,
            )
    """

    print("### COMPLETED GENE STITCH job ###")
    #print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
