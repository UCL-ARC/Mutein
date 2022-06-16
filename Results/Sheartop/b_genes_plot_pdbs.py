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
# not yet complete = ["TP53","PIK3CA","TP63"]

chosen_genes = ["NOTCH2"]
chosen_metric = "score"

def main(args):
            
    # current dir
    dir = os.path.dirname(os.path.realpath(__file__)) + "/Data/"
    
    for gene in chosen_genes:
        csv_file_bg = dir + gene.lower() + ".csv"                        
        df_bg = pd.read_csv(csv_file_bg)
        make_a_pdb_summary(gene,df_bg)    
    plt.show()
    


def make_a_pdb_summary(gene,df_in):        
    df = df_in.dropna(subset=["gene_no"])   
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    xax = "gene_no"
    yax = "pdb"
    hue = "ddg"    
    df["ddg"] = pd.to_numeric(df["ddg"])
    df["gene_no"] = pd.to_numeric(df["gene_no"])    
    df['pdb'] = df['pdb'].str.slice(0,-5) #assuming _rep10 on the end of them all :-)                        
    
    df = df.sort_values(by=[yax,"mut_from"], ascending=[False,True])           

    print(df)


    count = len(df.index)
    x_start = df["gene_no"].min()
    x_end = df["gene_no"].max()        
    dfG=df.groupby(['source','pdb','gene_no'])['ddg'].mean()
    dfG = dfG.reset_index()                
    fig.suptitle(gene + " PDB Gene Coverage\nddg <-1=stabilising >2.5=destabilising")                        
    vmin = -2.5  # ddg_df[hue].min() #they have defined >2.5 as destabilising
    vmax = -1 * vmin             
    dfG = dfG.sort_values(by=[yax], ascending=[False])                   
    g = ax.scatter(dfG[xax], dfG[yax], c=dfG[hue],cmap="Spectral",edgecolor="silver",alpha=0.35,linewidth=0.1,s=20,vmin=vmin,vmax=vmax)
    cb = plt.colorbar(g, extend="both")
    cb.set_label(hue)
    ax.set_xlabel(xax + "\n("+str(len(df.index))+")")
    ax.set_ylabel("")                 
    ax.xaxis.set_ticks(np.arange(0, x_end, 100)) 
    plt.xticks(rotation=90)   
    plt.tight_layout() 

        
    
if __name__ == "__main__":
    import sys
    globals()["main"](sys.argv)


