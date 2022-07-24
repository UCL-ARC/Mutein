"""
RSA 21/4/22
------------------------
Class to manage the submission to qsub on myriad

"""
import os
import subprocess
import Paths
import FileDf
from os.path import exists


class ScoringMetric:
    def __init__(self,gene_path,dataset,gene):
        self.dataset=dataset
        self.gene=gene
        self.path = gene_path
        pdb_file = self.path.outputs + "Coverage_all.csv"
        if exists(pdb_file):
            try:
                fdf = FileDf.FileDf(pdb_file,header=True)
                self.df = fdf.openDataFrame()        
                self.score_dict = {}                
                for i in range(len(self.df.index)):
                
                    src = self.df["source"][i]
                    pdb = self.df["pdb"][i]
                    mth = self.df["method"][i]
                    res = self.df["resolution"][i]
                    cov = self.df["coverage"][i]
                    #print(src,pdb,mth,res,cov)
                    metric = 0.8
                    if src == "SMHOM":
                        metric = 0.5
                    elif src == "AF":
                        metric = 0.01
                    elif src == "EXP" and "x-ray" in mth:
                        metric = 1 * 2/5/float(res)            
                    metric *= float(cov) * 1000
                    self.score_dict[pdb.lower()] = [round(metric,2),mth,res,cov] #score,method,resolution.coverage
            except:
                pass
            
    def cutPdb(self,pdb):
        if "_rep" in pdb:
            pos = pdb.find("_rep")
            pdb = pdb[:pos]
        return pdb        

    def getScore(self,pdb):
        if pdb.lower() in self.score_dict:
            return self.score_dict[pdb.lower()]
        else:
            return 0

