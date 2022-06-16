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

chosen_metric = "score"

def main(args):
            
    # current dir
    dir = os.path.dirname(os.path.realpath(__file__)) + "/Data/"            
    csv_file_summary = dir + "summary_stats_"+chosen_metric+".csv"
    df_summary = pd.read_csv(csv_file_summary)
    plot_summary_stats(df_summary,chosen_metric)                                
    plt.show()
    
def plot_summary_stats(df_summary,metric):     
    fig, [ax1,ax2,ax3] = plt.subplots(1, 3, figsize=(10, 5))

    #gene,
    # #all_bg_count,all_bg_mean,all_bg_median,all_var_count,all_var_mean,all_var_median,all_p_value,
    # mean_bg_count,mean_bg_mean,mean_bg_median,mean_var_count,mean_var_mean,mean_var_median,mean_p_value
    df_summary = df_summary.sort_values(by=["gene"],ascending=[False])
    
    fig.suptitle("Summary of statistical significance in ddg for gene set")                            

    xax,yax,hue="all_bg_median","all_var_median","all"
    ax1.scatter(df_summary[xax], df_summary[yax],label=hue,s=15,color="limegreen")
    xax,yax,hue=metric+"_bg_median",metric+"_var_median",metric
    ax1.scatter(df_summary[xax], df_summary[yax],label=hue,s=15,color="dodgerblue")
    ax1.set_xlabel("background median ddg")                         
    ax1.set_ylabel("variant median ddh")
    ax1.plot([0, 1], [0, 1],color="firebrick",linewidth=0.5)
    ax1.legend()

    xax,yax,hue="all_p_value","gene","all"
    ax2.scatter(df_summary[xax], df_summary[yax],label=hue,s=15,color="limegreen")
    xax,yax,hue=metric+"_p_value","gene",metric
    ax2.scatter(df_summary[xax], df_summary[yax],label=hue,s=15,color="dodgerblue")
    ax2.axvline(x=0.05,color="orangered",linewidth=0.5)
    ax2.set_xlabel("p-value")                         
    ax2.set_ylabel("gene")                         
    ax2.legend()
    
    binsnp = np.arange(0,1,0.025).tolist()    
    sns.histplot(df_summary["all_p_value"],alpha=0.5,bins=binsnp,label="all",ax=ax3,color="limegreen")
    sns.histplot(df_summary[metric+"_p_value"],alpha=0.6,bins=binsnp,label=metric,ax=ax3,color="dodgerblue")
    ax3.axvline(x=0.05,color="orangered",linewidth=0.5)
    
    ax3.legend()

    """
    xax,yax,hue,ax="all_bg_median","all_var_median","all_p_value",ax1
    g1 = ax.scatter(df_summary[xax], df_summary[yax], c=df_summary[hue],cmap="Spectral",edgecolor="silver",alpha=0.95,linewidth=0.1,s=20)
    cb1 = plt.colorbar(g1,ax=ax)
    cb1.set_label(hue)    
    ax.set_ylabel("")                         

    xax,yax,hue,ax=metric+"_bg_median",metric+"_var_median",metric+"_p_value",ax2
    g2 = ax.scatter(df_summary[xax], df_summary[yax], c=df_summary[hue],cmap="Spectral",edgecolor="silver",alpha=0.95,linewidth=0.1,s=20)
    cb2 = plt.colorbar(g2,ax=ax)
    cb2.set_label(hue)    
    ax.set_ylabel("")                         
    """
    
    
    plt.tight_layout() 
    
   
    
    


    
if __name__ == "__main__":
    import sys
    globals()["main"](sys.argv)


