

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
            return self.getPdbProgressReport(self.pdb)
        elif self.gene != "" and self.gene.lower() != "x":
            return self.getGeneProgressReport(self.gene)
        else:
            return self.createDatasetProgressReport(self.dataset)


    def createDatasetProgressReport(self,dataset):
        ds_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=dataset)                
        filename = ds_path.inputs + "genes_pdb_list.csv"        
        print("GENE      \t\tBACKGROUND  \t\tVARIANTS  ")
        print("----------\t\t----------\t\t----------")
        if exists(filename):
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for line in lines[1:]:
                    geneo = line.strip().split(",")[1]
                    patho = Paths.Paths(self.data_dir, self.pipe_dir, dataset=dataset, gene=geneo)
                    filenameA = patho.outputs + "ddg_bm_background.csv"
                    line_string = geneo+"\t\t"
                    existsfile, time = self.checkFile(filenameA)
                    if existsfile:
                        line_string += time + "\t\t"
                    else:
                        line_string += "----\t\t"
                    filenameB = patho.outputs + "ddg_variant_bm.csv"
                    existsfile, time = self.checkFile(filenameB)
                    if existsfile:
                        line_string += time + "\t\t"
                    else:
                        line_string += "----\t\t"
                    print(line_string)
        else:
            print("Gene has not been prepped")
            print("TODO: Submit pdb prepare")
        return ""

    def getGeneProgressReport(self,gene):
        gene_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene)                        
        filenameA = gene_path.outputs + "ddg_bm_background.csv"
        filenameB = gene_path.outputs + "ddg_variant_bm.csv"        
        existsfileA, timeA = self.checkFile(filenameA)
        existsfileB, timeB = self.checkFile(filenameB)
        print("SUMMARY BACKGROUND\t\tSUMMARY VARIANTS")
        print("------------------\t\t-----------------")
        exists_string = ""
        if existsfileA:
            exists_string += timeA + "\t\t"
        else:
            exists_string += "  ----  \t\t"
        if existsfileB:
            exists_string += timeB + "\t\t"
        else:
            exists_string += "  ----  \t\t"
        print(exists_string)        
        print("")        
        filename = gene_path.outputs + "pdb_tasklist.csv"
        if exists(filename):
            print("PDB\t\tPDB REPAIR\t\tBG SPLIT\t\tBG TASKS\t\tVAR SPLIT\t\tVAR TASKS")
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for ln in lines[1:]:
                    pdbo = ln.strip().split(",")[2]
                    patho = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset, gene=gene, pdb=pdbo)
                    
                    pdb_file = patho.inputs + pdbo.lower() + "_repx.pdb"
                    bg_split = patho.thruputs + "params_background.txt"
                    bg_results = patho.outputs + "ddg_background.csv"
                    var_split = patho.thruputs + "params_variants.txt"                     
                    var_results = patho.outputs + "ddg_buildmodel.csv"                     
                                                                      
                    existsfileA, timeResA = self.checkFile(pdb_file)
                    existsfileB, timeResB = self.checkFile(bg_split)
                    existsfileC, timeResC = self.checkFile(bg_results)
                    existsfileD, timeResD = self.checkFile(var_split)
                    existsfileE, timeResE = self.checkFile(var_results)

                    pdb_line = pdbo + "\t\t"                    
                    if existsfileA:
                        pdb_line += timeResA + "\t\t"
                    else:
                        pdb_line += " ---- \t\t"
                    if existsfileB:
                        pdb_line += timeResB + "\t\t"
                    else:
                        pdb_line += " ---- \t\t"
                    if existsfileC:
                        pdb_line += timeResC + "\t\t"
                    else:
                        pdb_line += self.completedPdbTaskFiles(gene,pdbo,False) +  "\t\t"
                    if existsfileD:
                        pdb_line += timeResD + "\t\t"
                    else:
                        pdb_line += " ---- \t\t"
                    if existsfileE:
                        pdb_line += timeResE + "\t\t"
                    else:
                        pdb_line += self.completedPdbTaskFiles(gene,pdbo,True) +  "\t\t"
                    print(pdb_line)
                        
        else:
            print("The pdbs have not been prepared")
            print("TODO: Submit pdb prepare")
        return ""
    
    def getPdbProgressReport(self,gene,pdb):
        pdb_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene,pdb=pdb)
        return ""
    
    def existsGenePdbFiles(self,gene):
        '''
        returns the number complete and the latest file stamp
        '''
        return ""
    
    def existsPdbFile(self,gene,pdb):
        '''
        returns if it is completed and the filestamp
        '''
        pdb_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene,pdb=pdb)
        return ""

    def existsGeneSplitFiles(self,gene):
        '''
        returns the number complete and the latest file stamp
        '''
        return ""
    
    def existsSplitFile(self,gene,pdb):
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
    
    def completedPdbTaskFiles(self,pdb,gene,isvariant):
        '''
        returns if it is completed and the filestamp
        '''
        pdb_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene,pdb=pdb)
        filenameo = pdb_path.thruputs + "params_background.txt"
        if isvariant:
            filenameo = pdb_path.thruputs + "params_variants.txt"
        existso, timeo = self.checkFile(filenameo)
        count = 0
        lasttime = ""                  
        msg = ""
        if existso:                            
            with open(filenameo, "r") as fr:
                lines = fr.readlines()
                numtasks = len(lines) - 1                                                                
                for i in range(1, len(lines)):
                    filenameoo = pdb_path.thruputs + "agg/" + str(i) + "_ddg_background.csv"
                    if isvariant:
                        filenameoo = pdb_path.thruputs + "vagg/" + str(i) + "_ddg_buildmodel.csv"
                    existst, timet = self.checkFile(filenameoo)
                    if existst:
                        count += 1
                        lasttime = timet                                  
        
        msg += str(count) + "/" + str(numtasks)        
        if count > 0:
            msg += "@" + str(lasttime)
        
        return msg

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

    

