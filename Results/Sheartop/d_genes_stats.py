"""
(RSA 13/6/22)
Sample script to show the contents of the csv files from Mutein for the
top genes from the shearwater dataset

"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

all_genes_list = ["NOTCH1","NOTCH2","FAT1","NOTCH3","ARID2","AJUBA","KMT2D","DICER1","NSD1","NOTCH4","CUL3","RB1"]
# not yet complete = ["TP53","PIK3CA","TP63"]

chosen_genes = ["NOTCH2"]
chosen_metric = "score"

def main(args):
            
    # current dir
    dir = os.path.dirname(os.path.realpath(__file__)) + "/Data/"
    
    for gene in chosen_genes:
        csv_file_bg = dir + gene.lower() + ".csv"
        csv_file_var = dir + gene.lower() + "_v.csv"                
        df_bg = pd.read_csv(csv_file_bg)
        df_var = pd.read_csv(csv_file_var)                
        make_ddg_stats(gene,"All", df_bg,df_var)

        csv_file_bg_metric = dir + gene.lower()+"_"+chosen_metric+".csv"
        csv_file_var_metric = dir + gene.lower() + "_v"+"_"+chosen_metric+".csv"        
        df_bg_metric = pd.read_csv(csv_file_bg_metric)
        df_var_metric = pd.read_csv(csv_file_var_metric)                                                                        
        make_ddg_stats(gene,chosen_metric,df_bg_metric,df_var_metric)
            
    plt.show()
    
def make_ddg_stats(gene,metric,df_bg_in,df_var_in):
    print("pdb_summary")
    fig, [ax1,ax2] = plt.subplots(1, 2, figsize=(8, 4))    
    df_bg = df_bg_in.dropna(subset=["gene_no"])
    df_var = df_var_in.dropna(subset=["gene_no"])          
    df_bg["ddg"] = pd.to_numeric(df_bg["ddg"])
    df_var["ddg"] = pd.to_numeric(df_var["ddg"])          
    df_bg = df_bg.sort_values(by=["ddg"], ascending=[True])                 
    df_var = df_var.sort_values(by=["ddg"], ascending=[True])                 
    
    

    # calculate stats
    count_bg = len(df_bg.index)
    count_var = len(df_var.index)        
    mean_bg = df_bg["ddg"].mean()
    mean_var = df_var["ddg"].mean()
    median_bg = df_bg["ddg"].median()
    median_var = df_var["ddg"].median()

    stats_bg = "mean=" + str(round(mean_bg,2))
    stats_bg += "\nmedian=" + str(round(median_bg,2))
    stats_var = "mean=" + str(round(mean_var,2))
    stats_var += "\nmedian=" + str(round(median_var,2))
    res = stats.mannwhitneyu(df_bg["ddg"],df_var["ddg"])
    U_p = res.pvalue
    stats_var += "\nMannWhitneyU p-val=" + str(round(U_p,4))
    



    fig.suptitle(gene + " (" + metric + ")" +" Background vs Variants")                        

    sns.histplot(data=df_bg, x="ddg", palette="tab20", ax=ax1, bins=50)
    ax1.set_ylabel("count")
    ax1.set_xlabel("Background ddg\ncount="+str(count_bg)+"\n"+stats_bg)

    sns.histplot(data=df_var, x="ddg", palette="tab20", ax=ax2, bins=50)
    ax2.set_ylabel("count")
    ax2.set_xlabel("Variant ddg\ncount="+str(count_var)+"\n"+stats_var)

    plt.tight_layout() 
             
    
if __name__ == "__main__":
    import sys
    globals()["main"](sys.argv)


