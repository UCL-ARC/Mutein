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
#all_genes_list = ["AJUBA"]
# not yet complete = ["TP53","PIK3CA","TP63"]

chosen_metric = "score"

def main(args):
            
    # current dir
    dir = os.path.dirname(os.path.realpath(__file__)) + "/Data/"

    dic_summary = {}
    dic_summary["gene"]=[]
    dic_summary["all_bg_count"]=[]
    dic_summary["all_bg_mean"]=[]
    dic_summary["all_bg_median"]=[]
    dic_summary["all_var_count"]=[]
    dic_summary["all_var_mean"]=[]
    dic_summary["all_var_median"]=[]
    dic_summary["all_p_value"]=[]
    dic_summary[chosen_metric+"_bg_count"]=[]
    dic_summary[chosen_metric+"_bg_mean"]=[]
    dic_summary[chosen_metric+"_bg_median"]=[]
    dic_summary[chosen_metric+"_var_count"]=[]
    dic_summary[chosen_metric+"_var_mean"]=[]
    dic_summary[chosen_metric+"_var_median"]=[]
    dic_summary[chosen_metric+"_p_value"]=[]
    
    for gene in all_genes_list:
        csv_file_bg = dir + gene.lower() + ".csv"
        csv_file_var = dir + gene.lower() + "_v.csv"                
        df_bg = pd.read_csv(csv_file_bg)
        df_var = pd.read_csv(csv_file_var)                
        csv_file_bg_metric = dir + gene.lower()+"_"+chosen_metric+".csv"
        csv_file_var_metric = dir + gene.lower() + "_v"+"_"+chosen_metric+".csv"        
        df_bg_metric = pd.read_csv(csv_file_bg_metric)
        df_var_metric = pd.read_csv(csv_file_var_metric)                                                                                
        dic_summary = add_ddg_stats(dic_summary,gene,df_bg,df_var,chosen_metric,df_bg_metric,df_var_metric)

    df = pd.DataFrame.from_dict(dic_summary)     
    df.to_csv(dir+"summary_stats_" + chosen_metric + ".csv",index=False)           
    
    print(df)
    
    
def add_ddg_stats(dic_summary,gene,df_bg_in,df_var_in,metric,df_bg_metric_in,df_var_metric_in):        
    df_bg = df_bg_in.dropna(subset=["gene_no"])
    df_var = df_var_in.dropna(subset=["gene_no"])          
    df_bg_metric = df_bg_metric_in.dropna(subset=["gene_no"])
    df_var_metric = df_var_metric_in.dropna(subset=["gene_no"])          
    df_bg["ddg"] = pd.to_numeric(df_bg["ddg"])
    df_var["ddg"] = pd.to_numeric(df_var["ddg"])          
    df_bg_metric["ddg"] = pd.to_numeric(df_bg_metric["ddg"])
    df_var_metric["ddg"] = pd.to_numeric(df_var_metric["ddg"])          
    df_bg = df_bg.sort_values(by=["ddg"], ascending=[True])                 
    df_var = df_var.sort_values(by=["ddg"], ascending=[True])                         
    df_bg_metric = df_bg_metric.sort_values(by=["ddg"], ascending=[True])                 
    df_var_metric = df_var_metric.sort_values(by=["ddg"], ascending=[True])                         
    # calculate stats
    count_bg = len(df_bg.index)
    count_var = len(df_var.index)        
    count_bg_metric = len(df_bg_metric.index)
    count_var_metric = len(df_var_metric.index)        
    mean_bg = df_bg["ddg"].mean()
    mean_var = df_var["ddg"].mean()
    median_bg = df_bg["ddg"].median()
    median_var = df_var["ddg"].median()
    mean_bg_metric = df_bg_metric["ddg"].mean()
    mean_var_metric = df_var_metric["ddg"].mean()
    median_bg_metric = df_bg_metric["ddg"].median()
    median_var_metric = df_var_metric["ddg"].median()    
    res = stats.mannwhitneyu(df_bg["ddg"],df_var["ddg"])
    res_metric = stats.mannwhitneyu(df_bg_metric["ddg"],df_var_metric["ddg"])
    U_p = res.pvalue
    U_p_metric = res_metric.pvalue
    
    dic_summary["gene"].append(gene)
    count_bg = len(df_bg.index)
    dic_summary["all_bg_count"].append(count_bg)
    dic_summary["all_bg_mean"].append(mean_bg)
    dic_summary["all_bg_median"].append(median_bg)
    dic_summary["all_var_count"].append(count_var)
    dic_summary["all_var_mean"].append(mean_var)
    dic_summary["all_var_median"].append(median_var)
    dic_summary["all_p_value"].append(U_p)
    dic_summary[metric+"_bg_count"].append(count_bg_metric)
    dic_summary[metric+"_bg_mean"].append(mean_bg_metric)
    dic_summary[metric+"_bg_median"].append(median_bg_metric)
    dic_summary[metric+"_var_count"].append(count_var_metric)
    dic_summary[metric+"_var_mean"].append(mean_var_metric)
    dic_summary[metric+"_var_median"].append(median_var_metric)
    dic_summary[metric+"_p_value"].append(U_p_metric)

    return dic_summary
    


    
if __name__ == "__main__":
    import sys
    globals()["main"](sys.argv)


