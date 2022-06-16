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

all_genes_list = ["NOTCH1","NOTCH2","FAT1","NOTCH3","ARID2","AJUBA","KMT2D","DICER1","NSD1","NOTCH4","CUL3","RB1"]
#all_genes_list = ["AJUBA"]
# not yet complete = ["TP53","PIK3CA","TP63"]

#all_genes_list = ["AJUBA"]

def main(args):
            
    # current dir
    dir = os.path.dirname(os.path.realpath(__file__)) + "/Data/"
    
    for gene in all_genes_list:
        csv_file_bg = dir + gene.lower() + ".csv"
        csv_file_var = dir + gene.lower() + "_v.csv"
        print(csv_file_bg,csv_file_var)        
        df_bg = pd.read_csv(csv_file_bg)
        df_var = pd.read_csv(csv_file_var)
        
        apply_metric_choice(dir,df_bg,gene,choice="score")
        apply_metric_choice(dir,df_var,gene+"_v",choice="score")

        #apply_metric_choice(dir,df_bg,gene,choice="newscore")
        #apply_metric_choice(dir,df_var,gene+"_v",choice="newscore")
                        
        apply_metric_choice(dir,df_bg,gene,choice="mean")
        apply_metric_choice(dir,df_var,gene+"_v",choice="mean")

                
def calculate_score(source,method,coverage,resolution,score):    
    # can mess about with the score here
    return score


def apply_metric_choice(dir,df_in,name,choice):    
    df = df_in
    df = df.dropna(subset=["gene_no"])
    df["ddg"] = pd.to_numeric(df["ddg"])
    df["gene_no"] = pd.to_numeric(df["gene_no"])
    print(df)   
    if choice == "newscore":  #apply own calculated metric    
        df['newscore'] = df.apply(lambda row: calculate_score(row['source'],row['method'],row['coverage'],row['resolution'],row['score']), axis=1)
         
    colls = ["pdb_chain","mut_from","mut_to","gene_no"]            
    if choice == "mean":
        df_metric= df.groupby(colls)['ddg'].mean()        
        df_metric = df_metric.reset_index()        
    else:
        metric = "score"
        colls.append("ddg")
        if choice == "newscore":  #apply own calculated metric                
            metric = "newscore"                            
        #df_metric = df.sort_values(metric).groupby(colls).tail(1)                        
        df_metric= df.groupby(colls)['score'].max()        
        df_metric = df_metric.reset_index()
        
    df_metric = df_metric.sort_values(by=["gene_no","mut_to"], ascending=[True,True])
    df_metric.to_csv(dir+name+"_"+choice+".csv",index=False)
    return df_metric




        
        
    
if __name__ == "__main__":
    import sys
    globals()["main"](sys.argv)


