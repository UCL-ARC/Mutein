"""
RSA 3/5/22
This performs required analysis for the ddg

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class Analysis:
    def __init__(self,df,pdb_gene):        
        self.df = df
        self.pdb_gene = pdb_gene
        self.count=len(self.df.index)
    
    def createDdgResidue(self,file_path,title,stabilising=-1, destabilising=2.5,xax="pdb_rid"):
        # And save something visual as a starting point for some analysis                
        self.df["ddg"] = pd.to_numeric(self.df["ddg"])
        df_stab =  self.df.query("ddg <= " + str(stabilising))        
        df_destab =  self.df.query("ddg >= " + str(destabilising))
        
        if len(df_stab.index)>0:
            str_ratio = "1:" + str(round(len(df_destab.index)/len(df_stab.index)))
        else:
            str_ratio = str(round(len(df_destab.index))) + ":0"
        
        str_ratio = str(round(len(df_stab.index))) + ":" + str(round(len(df_destab.index)))
                
        fig, (ax1,ax2,ax3) = plt.subplots(1, 3,figsize=(15,5))        
        fig.suptitle(self.pdb_gene + " " + title + " ("+str(self.count)+")\nddg <-1=stabilising >2.5=destabilising Ratio=" + str_ratio)                            
        yax = "mut_to"
        hue = "ddg"
        
        ###  first plt ######
        sns.histplot(data=self.df, x="ddg", palette="tab20", ax=ax1, bins=50)
        ax1.set_ylabel("")

        ###  second plt ######
        # Clip it to show stabilising and destabilising mutations
        capp = 5 #what do we want to cap out the histogram at to see it better                
        #df['a'][df['a'] >= maxVal] = maxVal this returns the annoying copy error
        np_clipped = np.clip(self.df['ddg'], a_max=capp, a_min=None)
        #sns.histplot(data=np_clipped, x="ddg", palette="tab20", ax=ax2, bins=[-5,-2.5,-1,0,1,2.5,5])
        sns.histplot(data=np_clipped, palette="tab20", ax=ax2, bins=[-5,-2.5,-1,0,1,2.5,5])
        ax2.set_ylabel("")
        ax2.set_xlabel("ddg max=" + str(capp))
            
        ###  third plt ######
        vmin = -2.5  # ddg_df[hue].min() #they have defined >2.5 as destabilising
        vmax = -1 * vmin
        ddg_df = self.df.sort_values(by=yax, ascending=False)
        g = ax3.scatter(
            ddg_df[xax],
            ddg_df[yax],
            c=ddg_df[hue],
            cmap="Spectral",
            edgecolor="silver",
            alpha=0.35,
            linewidth=0.1,
            s=20,
            vmin=vmin,
            vmax=vmax,
        )
        cb = plt.colorbar(g, extend="both")
        cb.set_label(hue)
        ax3.set_xlabel("residue no")
        ax3.set_ylabel("")
        # plt.legend(title=hue,bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.,shadow=False,fancybox=False)  # Put the legend out of the figure        
        plt.savefig(file_path)
        print("### Analysis:OutputDdgResidue to",file_path)

            