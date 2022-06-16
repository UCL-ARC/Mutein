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
chosen_metrics = ["score"]

def main(args):
            
    # current dir
    dir = os.path.dirname(os.path.realpath(__file__)) + "/Data/"
    
    for gene in chosen_genes:
        csv_file_bg = dir + gene.lower() + ".csv"
        csv_file_var = dir + gene.lower() + "_v.csv"                        
        df_bg = pd.read_csv(csv_file_bg)
        df_var = pd.read_csv(csv_file_var)        
        make_ddg_histogram(gene, "Background:All", df_bg)
        make_ddg_histogram(gene, "Variant:All", df_var)

        for chosen_metric in chosen_metrics:
            csv_file_bg_metric = dir + gene.lower()+"_"+chosen_metric+".csv"
            csv_file_var_metric = dir + gene.lower() + "_v"+"_"+chosen_metric+".csv"
            df_bg_metric = pd.read_csv(csv_file_bg_metric)
            df_var_metric = pd.read_csv(csv_file_var_metric)                                                                        
            make_ddg_histogram(gene, "Background:"+chosen_metric, df_bg_metric)
            make_ddg_histogram(gene, "Variant:"+chosen_metric, df_var_metric)

        
    plt.show()
    
def make_ddg_histogram(gene, name, df_in):    
    # And save something visual as a starting point for some analysis    
    df = df_in.dropna(subset=["gene_no"])
    df["ddg"] = pd.to_numeric(df["ddg"])        
    stab = -1
    destab = 2.5
    df_stab = df.query("ddg <= " + str(stab))
    df_destab = df.query("ddg >= " + str(destab))
    count = len(df.index)    
    str_ratio = str(round(len(df_stab.index))) + ":" + str(round(len(df_destab.index)))
    if len(df_stab.index) >0:
        str_ratio += " ≈ 1:" + str(round(len(df_destab.index)/len(df_stab.index)))

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(gene + " "+ name + " (" + str(count)+ ")\nddg <-1=stabilising >2.5=destabilising Ratio=" + str_ratio )
    yax = "mut_to"
    xax = "gene_no"
    hue = "ddg"

    ###  first plt ######
    sns.histplot(data=df, x="ddg", palette="tab20", ax=ax1, bins=50)
    ax1.set_ylabel("count")

    ###  second plt ####### Clip it to show stabilising and destabilising mutations
    capp = 5  
    np_clipped = np.clip(df["ddg"], a_max=capp, a_min=None)    
    binsnp = np.arange(-5, 5, 0.5).tolist()
    sns.histplot(data=np_clipped, palette="tab20", ax=ax2, bins=binsnp)#[-5, -2.5, -1, 0, 1, 2.5, 5])
    ax2.set_ylabel("count")
    ax2.set_xlabel("ddg capped at " + str(capp))

    ###  third plt ######
    vmin = -2.5  # ddg_df[hue].min() #they have defined >2.5 as destabilising
    vmax = -1 * vmin
    ddg_df = df.sort_values(by=yax, ascending=False)
    g = ax3.scatter(
        ddg_df[xax],
        ddg_df[yax],
        c=ddg_df[hue],
        cmap="Spectral",
        edgecolor="silver",
        alpha=0.65,
        linewidth=0.1,
        s=20,
        vmin=vmin,
        vmax=vmax,
    )
    cb = plt.colorbar(g, extend="both")
    cb.set_label(hue)
    ax3.set_xlabel(xax)
    ax3.set_ylabel("")       
    plt.tight_layout()                          


            
if __name__ == "__main__":
    import sys
    globals()["main"](sys.argv)


