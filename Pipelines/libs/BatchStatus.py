

"""
RSA 14/6/22
Analyses the status of a batch
-----------------------------


"""
import Paths
import os
import subprocess
import SubRunner
from os.path import exists
import pathlib
from datetime import datetime


class BatchStatus:
    def __init__(self, data_dir,pipe_dir,dataset, gene, pdb):
        self.dataset = dataset
        self.gene=gene
        self.pdb=pdb
        self.data_dir = data_dir
        self.pipe_dir = pipe_dir
                        
    def createReport(self):
        if self.pdb != "" and self.pdb.lower() != "x":
            return self.getPdbProgresReport(self.pdb)
        elif self.gene != "" and self.gene.lower() != "x":
            return self.getGeneProgresReport(self.gene)
        else:
            return self.getDatasetProgresReport(self.dataset)


    def createDatasetProgressReport(self,dataset):
        ds_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=dataset)                
        filename = ds_path.inputs + "genes_pdb_list.csv"        
        print("GENE    \tBACKGROUND\tVARIANTS")
        print("________\t___________\t________")
        if exists(filename):
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for line in lines[1:]:
                    geneo = line.strip().split(",")[1]
                    patho = Paths.Paths(self.data_dir, self.pipe_dir, dataset=dataset, gene=geneo)
                    filenameA = patho.outputs + "ddg_bm_background.csv"
                    line_string = "geneo\t"
                    existsfile, time = self.checkFile(filenameA)
                    if existsfile:
                        line_string += time + "\t"
                    else:
                        line_string += "----\t"
                    filenameB = patho.outputs + "ddg_variant_bm.csv"
                    existsfile, time = self.checkFile(filenameB)
                    if existsfile:
                        line_string += time + "\t"
                    else:
                        line_string += "----\t"
                    print(line_string)
        else:
            print("Gene has not been prepped")
            print("Run: Submit pdb prepare")
        return ""

    def getGeneProgressReport(self,gene):
        gene_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene)
        return ""
    
    def getPdbProgressReport(self,gene,pdb):
        pdb_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene,pdb=pdb)
        return ""
    
    def existsGenePdbFiles(self,gene):
        '''
        returns the number complete and the latest file stamp
        '''
        return ""
    
    def existsPdbFile(self,pdb):
        '''
        returns if it is completed and the filestamp
        '''
        return ""

    def existsGeneSplitFiles(self,gene):
        '''
        returns the number complete and the latest file stamp
        '''
        return ""
    
    def existsSplitFile(self,pdb):
        '''
        returns if it is completed and the filestamp
        '''
        return ""

    def completedGeneTaskFiles(self,gene,isvariant):
        '''
        returns the number complete and the latest file stamp
        AND those completed and those not
        '''
        return ""
    
    def completedPdbTaskFiles(self,pdb,isvariant):
        '''
        returns if it is completed and the filestamp
        '''
        return ""

    def completedGeneResultsFiles(self,gene,isvariant):
        '''
        returns the number complete and the latest file stamp
        AND those completed and those not
        '''
        return ""
    
    def completedPdbResultsFiles(self,pdb,isvariant):
        '''
        returns if it is completed and the filestamp
        '''
        return ""

    def checkFile(self,onefile):
        if exists(onefile):
            timeA = pathlib.Path(onefile).stat().st_mtime
            return True, datetime.fromtimestamp(timeA).strftime("%d%b%y-%H:%M")
        return False, ""

    

