"""
RSA 3/5/22
This performs required analysis for the ddg

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
pd.options.mode.chained_assignment = None  # default='warn'

class Analysis:
    def __init__(self, df_back, df_var, name):
        self.data_back = df_back
        self.data_var = df_var
        self.stabilising=-1
        self.destabilising=2.5
        self.name = name
        
    def createPdbSummary(self):
        #self.clearPlots()
        fig, ax = plt.subplots(1, 1, figsize=(10, 7))
        # And save something visual as a starting point for some analysis
        xax = "gene_no"
        yax = "pdb"
        hue = "ddg"
        df = self.data_back
        
        #mainpulate dataframe
        try:   
            df["ddg"] = pd.to_numeric(df["ddg"])
            df["gene_no"] = pd.to_numeric(df["gene_no"])
            df['pdb'] = df['pdb'].str.slice(0,-5) #assuming _rep10 on the end of them all :-)                
            df = df.dropna()     
            count = len(df.index)
            x_start = df["gene_no"].min()
            x_end = df["gene_no"].max()        
            df=df.groupby(['source','pdb','gene_no'])['ddg'].mean()
            df = df.reset_index()
            print(df)
            df['absddg'] = abs(df['ddg'])                
            fig.suptitle(self.name + " PDB Gene Coverage\nddg <-1=stabilising >2.5=destabilising")                
            ###  third plt ######
            vmin = -2.5  # ddg_df[hue].min() #they have defined >2.5 as destabilising
            vmax = -1 * vmin        
            #df = df.sort_values(by=[yax,"absddg"], ascending=[False,True])        
            df = df.sort_values(by=[yax,"gene_no"], ascending=[False,False])        
            g = ax.scatter(df[xax], df[yax], c=df[hue],cmap="Spectral",edgecolor="silver",alpha=0.35,linewidth=0.1,s=20,vmin=vmin,vmax=vmax)
            cb = plt.colorbar(g, extend="both")
            cb.set_label(hue)
            ax.set_xlabel(xax + "\n("+str(len(df.index))+")")
            ax.set_ylabel("")                 
            ax.xaxis.set_ticks(np.arange(0, x_end, 100)) 
            plt.xticks(rotation=90)   
            plt.tight_layout() 
        except:
            print("error in df")

    def createDdgBackRid(self):                
        xax="pdb_rid"                        
        df = self.data_back
        df["ddg"] = pd.to_numeric(df["ddg"])
        df["pdb_rid"] = pd.to_numeric(df["pdb_rid"])                
        count = int(len(df.index))
        fig, (ax1) = plt.subplots(1, 1, figsize=(10, 7))
        fig.suptitle(self.name + " Background Mutations - per residue\nddg <-1=stabilising >2.5=destabilising")
        yax = "mut_to"
        hue = "ddg"
        vmin = -2.5  # ddg_df[hue].min() #they have defined >2.5 as destabilising
        vmax = -1 * vmin
        df = df.sort_values(by=yax, ascending=False)
        print(df)
        #g = ax1.scatter(df["pdb_rid"],df["mut_to"],c=df["ddg"])
        g = ax1.scatter(df[xax],df[yax],c=df[hue],cmap="Spectral",edgecolor="silver",alpha=1,linewidth=0.1,s=20,vmin=vmin,vmax=vmax)
        cb = plt.colorbar(g, extend="both")
        cb.set_label(hue)        
        ax1.set_xlabel(xax + "\n(" + str(count)+")")
        plt.xticks(rotation=90)   
        ax1.set_ylabel("")
    
    def createDdgBackFromTo(self):                
        xax="mut_from"
        yax="mut_to"
        hue = "ddg"
        df = self.data_back
        df["ddg"] = pd.to_numeric(df["ddg"])        
        count = int(len(df.index))
        df=df.groupby(['mut_from','mut_to'])['ddg'].mean()
        df = df.reset_index()        
        fig, (ax1) = plt.subplots(1, 1, figsize=(10, 7))
        fig.suptitle(self.name + " Background mutations - From:To\nddg <-1=stabilising >2.5=destabilising")                
        #vmin,vmax = -2.5, 2.5  # ddg_df[hue].min() #they have defined >2.5 as destabilising
        vmin,vmax = -10, 10  # but for this plot we want more variation        
        df = df.sort_values(by=[yax,xax], ascending=[False,True])                
        g = ax1.scatter(df[xax],df[yax],c=df[hue],cmap="Spectral",edgecolor="silver",alpha=1,linewidth=0.1,s=20,vmin=vmin,vmax=vmax)
        cb = plt.colorbar(g, extend="both")
        cb.set_label(hue)
        plt.xticks(rotation=90)   
        ax1.set_xlabel(xax + "\n(" + str(count)+")")
        ax1.set_ylabel(yax)

    def histAllBackground(self):
        xax="ddg"                        
        df = self.data_back
        df["ddg"] = pd.to_numeric(df["ddg"])        
        count = int(len(df.index))
        fig, (ax1) = plt.subplots(1, 1, figsize=(10, 7))
        fig.suptitle(self.name + " Background mutations\nddg <-1=stabilising >2.5=destabilising")                
        vmin = -2.5  # ddg_df[hue].min() #they have defined >2.5 as destabilising
        vmax = -1 * vmin
        sns.histplot(data=df, x=xax, palette="tab20", ax=ax1, bins=50)
        ax1.set_ylabel("")                                 
        ax1.set_xlabel(xax + "\n(" + str(count)+")")
        ax1.set_ylabel("")
                                    
    def show(self):
        plt.show()

    def clearPlots(self):
        # clear all plots from memory        
        plt.figure().clear()
        plt.close()
        plt.cla()
        plt.clf()
