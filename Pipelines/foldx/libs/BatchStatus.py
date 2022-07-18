

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
            return self.getPdbProgressReport(self.gene,self.pdb)
        elif self.gene != "" and self.gene.lower() != "x":
            return self.getGeneProgressReport(self.gene)
        else:
            return self.createDatasetProgressReport(self.dataset)


    def getDirSize(self,start_path):
        total_size = 0        
        #for dirpath, dirnames, filenames in os.walk(start_path):
        #    for f in filenames:
        #        fp = os.path.join(dirpath, f)
        #        # skip if it is symbolic link
        #        if not os.path.islink(fp):
        #            total_size += os.path.getsize(fp)
        total_size += os.path.getsize(start_path)

        return str(total_size/1000000)+" MB"


    def createDatasetProgressReport(self,dataset):
        ds_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=dataset)                
        filename = ds_path.inputs + "genes_pdb_list.csv"        
        #print("GENE      \t\tSIZE       \t\tPDBS   \t\tBACKGROUND  \t\tVARIANTS  ")
        #print("----------\t\t-----------\t\t-------\t\t------------\t\t----------")
        print("GENE      \t\tPDBS   \t\tBACKGROUND  \t\tVARIANTS  ")
        print("-----------\t\t-------\t\t------------\t\t----------")
        if exists(filename):
            with open(filename, "r") as fr:
                lines = fr.readlines()
                genes = []
                for line in lines[1:]:
                    geneo = line.strip().split(",")[1]
                    genes.append(geneo)
                genes.sort()

                for geneo in genes:
                    geneo = line.strip().split(",")[1]
                    ####################
                    line_string = geneo+"  \t\t"
                    patho = Paths.Paths(self.data_dir, self.pipe_dir, dataset=dataset, gene=geneo)
                    ## dir size                    
                    #try:
                    #    size = self.getDirSize(patho.gene_inputs)
                    #    line_string += size + "\t\t      "
                    #except:
                    #    size = "xyz MB"
                    #    line_string += size + "\t\t      "
                    ###################
                    filenameC = patho.outputs + "pdb_tasklist.csv"                                        
                    existsfile, time = self.checkFile(filenameC)
                    if existsfile:
                        num_pdbs = self.getGeneNumPdbs(geneo)
                        line_string += num_pdbs + "\t\t"
                    else:                        
                        line_string += "----\t\t"                        
                    ####################
                    filenameA = patho.outputs + "ddg_background.csv"
                    filenameA_ds = ds_path.outputs + f"{geneo}_ddg_background.csv"
                    
                    existsfile, time = self.checkFile(filenameA_ds)
                    if existsfile:
                        line_string += "DS:"+time + "\t\t"
                    else:
                        existsfile, time = self.checkFile(filenameA)
                        if existsfile:
                            line_string += "G:"+time + "\t\t"
                    
                    if not existsfile:
                        num_done,num_tsks = self.getGeneNumTasks(geneo,False)
                        if num_tsks == 0:                        
                            line_string += "----\t\t"
                        elif num_tsks == num_done:
                            line_string += f"{num_tsks}\t\t"      
                        else:
                            line_string += f"{num_done}/{num_tsks}\t\t"      
                    #######################
                    filenameB = patho.outputs + "ddg_variants.csv"
                    filenameB_ds = ds_path.outputs + f"{geneo}_ddg_variants.csv"    
                    existsfile, time = self.checkFile(filenameB_ds)
                    if existsfile:
                        line_string += "DS:"+time + "\t\t"
                    else:
                        existsfile, time = self.checkFile(filenameB)
                        if existsfile:
                            line_string += "Ge:"+time + "\t\t"
                    
                    if not existsfile:                    
                        num_done,num_tsks = self.getGeneNumTasks(geneo,True)
                        if num_tsks == 0:                        
                            line_string += "----\t\t"
                        elif num_tsks == num_done:
                            line_string += f"{num_tsks}\t\t"      
                        else:
                            line_string += f"{num_done}/{num_tsks}\t\t"      
                    print(line_string)
        else:
            print("Gene has not been prepped")
            print("TODO: Submit pdb prepare")
        return ""

    def getGenes(self):
        genes = []
        ds_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset)                
        filename = ds_path.inputs + "genes_pdb_list.csv"                
        if exists(filename):
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for line in lines[1:]:
                    geneo = line.strip().split(",")[1]
                    genes.append(geneo)                    
        return genes
    
    def getGeneProgressReport(self,gene):
        gene_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene)                        
        filenameA = gene_path.outputs + "ddg_background.csv"
        filenameB = gene_path.outputs + "ddg_variant.csv"        
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
                    if pdbo[0] != "#":
                        patho = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset, gene=gene, pdb=pdbo)
                        
                        pdb_file = patho.inputs + pdbo.lower() + "_repx.pdb"
                        bg_split = patho.inputs + "params_background.txt"
                        bg_results = patho.outputs + "ddg_background.csv"
                        var_split = patho.inputs + "params_variants.txt"                     
                        var_results = patho.outputs + "ddg_variants.csv"                     
                                                                        
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
                            numC = self.getPdbNumTasks(gene,pdbo,False)
                            #pdb_line += numC+"\t\t"
                            pdb_line += timeResC+"\t\t"
                        else:
                            msgC,tmC = self.completedPdbTaskFiles(gene,pdbo,False)
                            pdb_line += msgC + "\t\t"
                            #pdb_line += tmC + "\t\t"
                        if existsfileD:
                            pdb_line += timeResD + "\t\t"
                        else:
                            pdb_line += " ---- \t\t"                    
                        if existsfileE:
                            numE = self.getPdbNumTasks(gene,pdbo,True)
                            #pdb_line += numE+"\t\t"
                            pdb_line += timeResE+"\t\t"
                        else:
                            msgE,tmE = self.completedPdbTaskFiles(gene,pdbo,True)
                            pdb_line += msgE + "\t\t"
                            #pdb_line += tmE + "\t\t"
                        print(pdb_line)
                        
        else:
            print("The pdbs have not been prepared")
            print("TODO: Submit pdb prepare")
        return ""
    
    def getGenePdbs(self,gene):
        gene_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene)
        filename = gene_path.outputs + "pdb_tasklist.csv"
        pdbs = []
        if exists(filename):            
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for ln in lines[1:]:
                    pdbo = ln.strip().split(",")[2]
                    if pdbo[0] != "#":
                        pdbs.append(pdbo)
        return pdbs
    
    def getGeneNumPdbs(self,gene):
        gene_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene)
        filename = gene_path.outputs + "pdb_tasklist.csv"
        num_pdbs = 0
        num_done = 0
        if exists(filename):            
            with open(filename, "r") as fr:
                lines = fr.readlines()
                if len(lines) > 0:                    
                    for ln in lines[1:]:
                        pdbo = ln.strip().split(",")[2]
                        if pdbo[0] != "#":
                            num_pdbs += 1
                            patho = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset, gene=gene, pdb=pdbo)                        
                            pdb_file = patho.inputs + pdbo.lower() + "_repx.pdb"
                            existsfile, timeRes = self.checkFile(pdb_file)
                            if existsfile:
                                num_done += 1
                                        
            if num_done == num_pdbs:
                pdb_summary = str(num_done)
            else:
                pdb_summary = str(num_done) + "/" + str(num_pdbs)
            return str(pdb_summary)
        else:
            return "-"                

    def getPdbNumTasksComplete(self,gene,pdb,isvariant):
        pdb_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene,pdb=pdb)
        filenameo = pdb_path.inputs + "params_background.txt"
        if isvariant:
            filenameo = pdb_path.inputs + "params_variants.txt"
        existso, timeo = self.checkFile(filenameo)        
        nums = 0
        if existso:                            
            with open(filenameo, "r") as fr:
                lines = fr.readlines()
                numtasks = len(lines)-1
                for n in range(numtasks):
                    no = n+1
                    filenamet = pdb_path.thruputs + f"agg/{no}_ddg_background.csv"
                    if isvariant:
                        filenamet = pdb_path.thruputs + f"vagg/{no}_ddg_variants.csv"
                    existst, timeo = self.checkFile(filenamet)      
                    if existst:
                        nums+=1  
        return nums

    def getGeneNumTasks(self,gene,isvariant):
        gene_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene)
        filename = gene_path.outputs + "pdb_tasklist.csv"
        pdbs = []
        num_tasks = 0
        num_done = 0
        if exists(filename):            
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for ln in lines[1:]:
                    pdbo = ln.strip().split(",")[2]
                    if pdbo[0] != "#":
                        num_tasks += int(self.getPdbNumTasks(gene,pdbo,isvariant))
                        num_done += self.getPdbNumTasksComplete(gene,pdbo,isvariant)
                    
        return num_done,num_tasks

    def getPdbNumTasks(self,gene,pdb,isvariant):
        pdb_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene,pdb=pdb)
        filenameo = pdb_path.inputs + "params_background.txt"
        if isvariant:
            filenameo = pdb_path.inputs + "params_variants.txt"
        existso, timeo = self.checkFile(filenameo)        
        numtasks = "0"
        if existso:                            
            with open(filenameo, "r") as fr:
                lines = fr.readlines()
                numtasks = str(len(lines)-1)
        return numtasks

    def getPdbProgressReport(self,gene,pdb):
        pdb_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene,pdb=pdb)                        
        filenameA = pdb_path.inputs + pdb + "_repx.pdb"
        filenameB = pdb_path.outputs + "ddg_background.csv"
        filenameC = pdb_path.outputs + "ddg_variants.csv"                        
        filenameD = pdb_path.inputs + "params_background.txt"
        filenameE = pdb_path.inputs + "params_variants.txt"
        
        existsfileA, timeA = self.checkFile(filenameA)
        existsfileB, timeB = self.checkFile(filenameB)
        existsfileC, timeC = self.checkFile(filenameC)
        existsfileD, timeD = self.checkFile(filenameD)
        existsfileE, timeE = self.checkFile(filenameE)
                
        print("PDB REPAIR\t\tPDB BACKGROUND\t\tPDB VARIANTS\t\tBG SPLIT\t\tVAR SPLIT")
        print("-----------\t\t--------------\t\t------------\t\t---------\t\t---------")
        exists_string = ""
        if existsfileA:
            exists_string += timeA + "\t\t"
        else:
            exists_string += "  ----  \t\t"
        if existsfileB:
            exists_string += timeB + "\t\t"
        else:
            exists_string += "  ----  \t\t"
        if existsfileC:
            exists_string += timeC + "\t\t"
        else:
            exists_string += "  ----  \t\t"
        if existsfileD:
            exists_string += timeD + "\t\t"
        else:
            exists_string += "  ----  \t\t"
        if existsfileE:
            exists_string += timeE + "\t\t"
        else:
            exists_string += "  ----  \t\t"
        print(exists_string)        
        print("")     

        if not existsfileA:
            print("TODO Submit repair")
        
        if not existsfileD or not existsfileE:
            print("TODO Submit splits prepare")
                
        if filenameD:
            count = 0
            print("------------- Background tasks -------------")
            with open(filenameD, "r") as fr:
                lines = fr.readlines()                                
                for i in range(1, len(lines)):
                    filenameo = pdb_path.thruputs + "agg/" + str(i) + "_ddg_background.csv"
                    existsfile, time = self.checkFile(filenameo)
                    if existsfile:
                        count += 1
                        print("Task", str(i), "at", time)
                    else:
                        print("Task", str(i), "----")
            print("------------------")
            print("Completed", count, "out of", len(lines) - 1)
            print("------------------")
            print("")
        
        if filenameE:
            count = 0
            print("------------- Variant tasks -------------")
            with open(filenameE, "r") as fr:
                lines = fr.readlines()                                
                for i in range(1, len(lines)):
                    filenameo = pdb_path.thruputs + "vagg/" + str(i) + "_ddg_variants.csv"
                    existsfile, time = self.checkFile(filenameo)
                    if existsfile:
                        count += 1
                        print("Task", str(i), "at", time)
                    else:
                        print("Task", str(i), "----")
            print("------------------")
            print("Completed", count, "out of", len(lines) - 1)
            print("------------------")

        

    def createUntasksForPdb(self,gene,pdb,write_genes=False,gene_writer=None,gene_var_writer=None,gene_num=0):
        path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene,pdb=pdb)        
        print("RECREATING TASK FILE with missing tasks\n")
        print("Check results files for pdb",pdb)
        #print(path.outputs)
        filenameA = path.outputs + "ddg_background.csv"
        filenameB = path.outputs + "ddg_variants.csv"       
        self.checkFile(filenameA) 
        self.checkFile(filenameB) 
        
        # Check the background
        filenameP = path.inputs + "params_background.txt"
        filename_incomplete = path.inputs + "params_background_incomplete.txt"
        count = 0
        print("\nChecking the background tasks")
        if exists(filenameP):
            with open(filename_incomplete, "w") as fw:
                with open(filenameP, "r") as fr:
                    lines = fr.readlines()                    
                    fw.write((lines[0]).strip() + "\n")
                    if write_genes and gene_num == 0:
                        gene_writer.write((lines[0]).strip() + "\n")
                    
                    for i in range(1, len(lines)):
                        filenameo = (
                            path.thruputs + "agg/" + str(i) + "_ddg_background.csv"
                        )
                        existsfile, time = self.checkFile(filenameo)
                        if existsfile:
                            count += 1                    
                        else:                            
                            fw.write((lines[i]).strip() + "\n")
                            if write_genes:
                                gene_writer.write((lines[i]).strip() + "\n")

            print("Completed", count, "out of", len(lines) - 1)

        else:
            print("Missing parameters file, the data needs preparation")

        # Check the background
        filenameP = path.inputs + "params_variants.txt"
        filename_incomplete = path.inputs + "params_variants_incomplete.txt"
        print("\nChecking the variant tasks")
        count = 0
        if exists(filenameP):
            with open(filename_incomplete, "w") as fw:
                with open(filenameP, "r") as fr:
                    lines = fr.readlines()
                    fw.write((lines[0]).strip() + "\n")
                    if write_genes and gene_num == 0:
                        gene_var_writer.write((lines[0]).strip() + "\n")                    
                    for i in range(1, len(lines)):
                        filenameo = (
                            path.thruputs + "vagg/" + str(i) + "_ddg_variants.csv"
                        )
                        existsfile, time = self.checkFile(filenameo)
                        if existsfile:
                            count += 1                            
                        else:                            
                            fw.write((lines[i]).strip() + "\n")
                            if write_genes:
                                gene_var_writer.write((lines[i]).strip() + "\n")                    
            print("Completed", count, "out of", len(lines) - 1)

        else:
            print(
                "Missing variants file, the data needs preparation, or there are none"
            )
    
                
    def makeMissingPdbs(self,gene):
        gene_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene)
        filename = gene_path.outputs + "pdb_tasklist.csv"        
        filename_incomplete = gene_path.outputs + "pdb_tasklist_incomplete.csv"              
        if exists(filename):
            with open(filename_incomplete, "w") as fw:
                with open(filename, "r") as fr:
                    lines = fr.readlines()
                    for ln in lines:
                        pdb = ln.strip().split(",")[2]
                        if pdb[0]!= "#":
                            if not self.existsPdbFile(gene,pdb):
                                fw.write((ln).strip() + "\n")

        
    def existsPdbFile(self,gene,pdb):
        '''
        returns if it is completed and the filestamp
        '''
        pdb_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene,pdb=pdb)
        filenameo = pdb_path.inputs + pdb.lower() + "_repx.pdb"
        exists,dt = self.checkFile(filenameo)
        return exists        
       
    def completedPdbTaskFiles(self,gene,pdb,isvariant):
        '''
        returns if it is completed and the filestamp
        '''
        pdb_path = Paths.Paths(self.data_dir, self.pipe_dir, dataset=self.dataset,gene=gene,pdb=pdb)
        filenameo = pdb_path.inputs + "params_background.txt"
        if isvariant:
            filenameo = pdb_path.inputs + "params_variants.txt"
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
                        filenameoo = pdb_path.thruputs + "vagg/" + str(i) + "_ddg_variants.csv"
                    existst, timet = self.checkFile(filenameoo)
                    if existst:
                        count += 1
                        lasttime = timet                                  
        
            if count == numtasks:
                msg += str(numtasks)        
            else:
                msg += str(count) + "/" + str(numtasks)        
                    
        return msg,lasttime


    def checkFile(self,onefile):
        if exists(onefile):
            timeA = pathlib.Path(onefile).stat().st_mtime
            return True, datetime.fromtimestamp(timeA).strftime("%d%b%y-%H:%M")
        return False, ""

    

