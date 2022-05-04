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



def run_pipeline(args):
    print("### Foldx aggregate pos scan ###")
    print(args)
    ##############################################
    argus = Arguments.Arguments(args)
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    gene_path = Paths.Paths("geneprot",dataset=dataset,gene=gene)    
    work_path = gene_path.gene_outputs + "agg/"
    argus.params["work_path"] = work_path
    gene_path.goto_job_dir(argus.arg("work_path"), args, argus.params, "_inputs01")
    ############################################
    params_file = gene_path.gene_outputs + "pdb_coverage.csv"
    fdfp = FileDf.FileDf(params_file,cols=["gene","accession","pdb","chain","method","resolution","coverage"])
    pm_df = fdfp.openDataFrame()
    all_df_vars = []
    all_df_backs = []
    all_df_vars_AF = []
    all_df_backs_AF = []
    for i in range(len(pm_df.index)):
        r = pm_df["pdb"][i].lower()
        # the file has already been turned into a dataframe called posscan_df.csv        
        pdb_path = Paths.Paths("pdb",dataset=dataset,gene=gene,pdb=r)    
        file_var = pdb_path.pdb_outputs + "ddg_variants.csv"
        file_back = pdb_path.pdb_outputs + "ddg_background.csv"
        if exists(file_var):
            fdf = FileDf.FileDf(file_var)
            all_df_vars.append(fdf.openDataFrame())        
            if r[:2] == "AF":
                all_df_vars_AF.append(fdf.openDataFrame())    
            else:
                all_df_vars.append(fdf.openDataFrame())        
        if exists(file_back):
            fdf = FileDf.FileDf(file_back)
            if r[:2] == "AF":
                all_df_backs_AF.append(fdf.openDataFrame())        
            else:
                all_df_backs.append(fdf.openDataFrame())        
            
    # NON Alpha-Fold structures
    outputs = []
    outputs.append([all_df_vars,"ddg_variants.csv","ddg_variants.png","variants","gene_no"])
    outputs.append([all_df_backs,"ddg_background.csv","ddg_background.png","background","pdb_rid"])
    outputs.append([all_df_vars_AF,"ddg_variants_AF.csv","ddg_variants_AF.png","variants AF","gene_no"])
    outputs.append([all_df_backs_AF,"ddg_background_AF.csv","ddg_background_AF.png","background AF","pdb_rid"])

    for dfs,csv_file, png_file, title,xax in outputs:
        if len(dfs) > 0:
            ddg_df = pd.concat(dfs,ignore_index=True)            
            df_file = gene_path.gene_outputs + csv_file
            ddg_df.to_csv(df_file, index=False)            
            plot_file = gene_path.gene_outputs + png_file
            anav = Analysis.Analysis(ddg_df,gene)
            anav.createDdgResidue(plot_file,title,xax=xax)
    

    #ddg_df_var = pd.concat(all_df_vars,ignore_index=True)        
    #ddg_df_back = pd.concat(all_df_backs,ignore_index=True)        
    #df_file_var = gene_path.gene_outputs + "ddg_variants.csv"
    #df_file_back = gene_path.gene_outputs + "ddg_background.csv"    
    #ddg_df_var.to_csv(df_file_var, index=False)        
    #ddg_df_back.to_csv(df_file_back, index=False)        
    #plot_file_var = gene_path.gene_outputs + "ddg_variants.png"
    #plot_file_back = gene_path.gene_outputs + "ddg_background.png"
    #anav = Analysis.Analysis(ddg_df_var,"")
    #anav.createDdgResidue(plot_file_var,"variant")
    #anab = Analysis.Analysis(ddg_df_var,"")
    #anab.createDdgResidue(plot_file_back,"background")
    
    print("### COMPLETED GENE STITCH job ###")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
