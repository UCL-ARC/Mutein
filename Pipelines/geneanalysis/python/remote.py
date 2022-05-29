"""
-----------------------------
RSA 24/05/22
-----------------------------

This simply views a file
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
import os
import pathlib
import pwd
import pandas as pd
from shutil import copyfile
from os.path import exists
import subprocess
from datetime import datetime

# import from the shared library in Mutein/Pipelines/shared/lib
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)

import Arguments
import Paths
from os import listdir
from os.path import isfile, join
import os


def checkResult(onefile):
    if exists(onefile):
        timeA = pathlib.Path(onefile).stat().st_mtime
        return True,datetime.fromtimestamp(timeA).strftime('%d-%m-%y-%H:%M')
    return False,""

def checkResults(ddg,bm,ps):
    if exists(ddg):
        timeA = pathlib.Path(ddg).stat().st_mtime            
        print("DDG Background file was created at",datetime.fromtimestamp(timeA).strftime('%d-%m-%y-%H:%M'))
    else:
        print("!!!DDG Background file does not exist")
    if exists(bm):
        timeB = pathlib.Path(bm).stat().st_mtime
        print("Variant buildmodel file was created at",datetime.fromtimestamp(timeB).strftime('%d-%m-%y-%H:%M'))
    else:
        print("!!!Variant buildmodel file does not exist")
    if exists(ps):
        timeC = pathlib.Path(ps).stat().st_mtime
        print("Variant posscan file was created at",datetime.fromtimestamp(timeC).strftime('%d-%m-%y-%H:%M'))
    else:
        print("!!!Variant posscan file does not exist")



def run_pipeline(args):    
    now = datetime.now()
    current_time = now.strftime("%d-%m-%y@%H.%H.%S")
    print("Mutein remote script:",current_time)
    
    ##############################################        
    mode = args[1]
    pattern = args[2]
    WorkDir = args[3]
    DataDir = args[4]
    InstallDir = args[5]
    PipelineDir = args[6]

    print("Mode=",mode)
    print("Pattern=",pattern)
    print("")
    if mode == "GENES":
        dataset_gene_pdb=pattern.split(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset)
        filename = path.inputs + "genes_pdb_list.csv"
        print("\nCheck genes list\n")      
        if exists(filename):
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for line in lines[1:]:                    
                    geneo =line.strip().split(",")[1]
                    patho = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=geneo)
                    filenameA = patho.outputs + "ddg_background.csv"
                    existsfile,time = checkResult(filenameA)
                    if existsfile:
                        print(geneo,time)
                    else:
                        print(geneo,"---")
        else:
            print("The dataset has not been prepared - no genes list",filename)
    elif mode == "PDBS":
        dataset_gene_pdb=pattern.split(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene)
        
        print("Check results files for gene")
        print(path.outputs)
        filenameA = path.outputs + "ddg_background.csv"
        filenameB = path.outputs + "ddg_variant_bm.csv"
        filenameC = path.outputs + "ddg_variant_ps.csv"
        checkResults(filenameA,filenameB,filenameC)

        print("\nCheck pdb list\n")        
        filename = path.outputs + "pdb_tasklist.csv"
        if exists(filename):
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for ln in lines[1:]:
                    pdbo =ln.strip().split(",")[2]
                    patho = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene,pdb=pdbo)
                    filenameA = patho.outputs + "ddg_background.csv"
                    filenamePdb = patho.inputs + pdb+"_rep10.pdb"
                    existsfile,time = checkResult(filenameA)
                    existsfilePdb,timePdb = checkResult(filenamePdb)
                    if existsfile and existsfilePdb:
                        print(pdbo,"\tPdb ready at", timePdb, "\tSplit ready at", time)
                    elif existsfile:
                        print(pdbo,"\tPdb not ready ---\tSplit ready at", time)
                    elif existsfilePdb:
                        print(pdbo,"\tPdb ready at", timePdb, "\tSplit not ready ---")
                    else:
                        print(pdbo,"--- --- --- ---")
                    
        else:
            print("The pdbs have not been prepared - no pdb list",filename)
    elif mode == "PDB":
        dataset_gene_pdb=pattern.split(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene,pdb=pdb)
        print("Check results files for pdb")
        print(path.outputs)
        filenameA = path.outputs + "ddg_background.csv"
        filenameB = path.outputs + "ddg_buildmodel.csv"
        filenameC = path.outputs + "ddg_posscan.csv"
        checkResults(filenameA,filenameB,filenameC)

        # check the pdb
        filenamePdb = path.inputs + pdb + "_rep10.pdb"        
        existsPdb,time = checkResult(filenamePdb)
        print("\nChecking the pdb repair")
        if existsPdb:
            print("...PDB 10 repair at",time)       
        else:
            print("...PDB 10 repair has not been done")       
        # Check the background        
        filenameP = path.thruputs + "params_background.txt"        
        existsP,time = checkResult(filenameP)
        count = 0
        print("\nChecking the background tasks")
        if existsP:     
            print("...Params background at",time)       
            with open(filenameP, "r") as fr:
                lines = fr.readlines()                
                print("\nThe pdb has been split into tasks=",len(lines)-1)
                print("...Any tasks that have completed are below\n")                
                for i in range(1,len(lines)):
                    filenameo = path.thruputs + "agg/" + str(i) + "_ddg_background.csv"
                    existsfile,time = checkResult(filenameo)
                    if existsfile:
                        count += 1
                        print("Task",str(i),"at",time)
                    else:
                        print("Task",str(i),"----")            
            print("Completed",count,"out of",len(lines)-1)
                    
        else:
            print("Missing parameters file, the data needs preparation")
        
        # Check the background
        filenameP = path.thruputs + "params_variants.txt"        
        print("\nChecking the variant tasks")
        count = 0
        if exists(filenameP):            
            with open(filenameP, "r") as fr:
                lines = fr.readlines()                    
                print("The variants have been split into tasks=",len(lines)-1)
                print("...Any tasks that have completed are below\n")                
                for i in range(1,len(lines)):
                    filenameo = path.thruputs + "vagg/" + str(i) + "_ddg_buildmodel.csv"
                    existsfile,time = checkResult(filenameo)
                    if existsfile:
                        count += 1
                        print("Task",str(i),"at",time)
                    else:
                        print("Task",str(i),"----")                        
            print("Completed",count,"out of",len(lines)-1)                    
        else:
            print("Missing variants file, the data needs preparation, or there are none")
    elif mode == "PDBINCOMPLETE":
        dataset_gene_pdb=pattern.split(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene,pdb=pdb)
        print("RECREATING TASK FILE with missing tasks\n")
        print("Check results files for pdb")
        print(path.outputs)
        filenameA = path.outputs + "ddg_background.csv"
        filenameB = path.outputs + "ddg_buildmodel.csv"
        filenameC = path.outputs + "ddg_posscan.csv"
        checkResults(filenameA,filenameB,filenameC)

        # Check the background
        filenameP = path.thruputs + "params_background.txt"
        filename_incomplete = path.thruputs + "params_background_incomplete.txt"
        count = 0
        print("\nChecking the background tasks")
        if exists(filenameP):
            with open(filename_incomplete, "w") as fw:
                with open(filenameP, "r") as fr:
                    lines = fr.readlines()
                    fw.write((lines[0]).strip()+"\n")
                    print("The pdb has been split into tasks=",len(lines)-1)
                    print("...Any tasks that have completed are below\n")                
                    for i in range(1,len(lines)):
                        filenameo = path.thruputs + "agg/" + str(i) + "_ddg_background.csv"
                        existsfile,time = checkResult(filenameo)
                        if existsfile:
                            count += 1
                            print("Task",str(i),"at",time)
                        else:
                            print("Task",str(i),"----")
                            fw.write((lines[i]).strip()+"\n")
            print("Completed",count,"out of",len(lines)-1)
                    
        else:
            print("Missing parameters file, the data needs preparation")
        
        # Check the background
        filenameP = path.thruputs + "params_variants.txt"
        filename_incomplete = path.thruputs + "params_variants_incomplete.txt"
        print("\nChecking the variant tasks")
        count = 0
        if exists(filenameP):
            with open(filename_incomplete, "w") as fw:
                with open(filenameP, "r") as fr:
                    lines = fr.readlines()
                    fw.write((lines[0]).strip()+"\n")
                    print("The variants have been split into tasks=",len(lines)-1)
                    print("...Any tasks that have completed are below\n")                
                    for i in range(1,len(lines)):
                        filenameo = path.thruputs + "vagg/" + str(i) + "_ddg_buildmodel.csv"
                        existsfile,time = checkResult(filenameo)
                        if existsfile:
                            count += 1
                            print("Task",str(i),"at",time)
                        else:
                            print("Task",str(i),"----")
                            fw.write((lines[i]).strip()+"\n")
            print("Completed",count,"out of",len(lines)-1)
                    
        else:
            print("Missing variants file, the data needs preparation, or there are none")
    elif mode == "PDB_BACK":
        dataset_gene_pdb=pattern.split(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene,pdb=pdb)        
        filename = path.outputs + "ddg_background.csv"        
        print(filename)
        mexists, time = checkResult(filename)
        if mexists:
            print("DATAFRAME_START")
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for line in lines:                
                    print(line.strip())
            print("DATAFRAME_END")
    elif mode == "PDB_BM":
        dataset_gene_pdb=pattern.split(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene,pdb=pdb)                
        filename = path.outputs + "ddg_buildmodel.csv"        
        if pdb == "":
            filename = path.outputs + "ddg_variant_bm.csv"
        print(filename)
        mexists, time = checkResult(filename)
        if mexists:
            print("DATAFRAME_START")
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for line in lines:                
                    print(line.strip())
            print("DATAFRAME_END")
    elif mode == "PDB_PS":
        dataset_gene_pdb=pattern.split(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene,pdb=pdb)                
        filename = path.outputs + "ddg_posscan.csv"        
        if pdb == "":
            filename = path.outputs + "ddg_variant_ps.csv"
            print(filename)
        mexists, time = checkResult(filename)
        if mexists:
            print("DATAFRAME_START")
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for line in lines:                
                    print(line.strip())
            print("DATAFRAME_END")
    elif mode == "COVERAGE":
        dataset_gene_pdb=pattern.split(":")
        dataset,gene,pdb = dataset_gene_pdb[0],dataset_gene_pdb[1],dataset_gene_pdb[2]
        path = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene)                
        filename = path.outputs + "Coverage_all.csv"
        mexists, time = checkResult(filename)
        if mexists:
            print("DATAFRAME_START")
            with open(filename, "r") as fr:
                lines = fr.readlines()
                for line in lines:                
                    print(line.strip())
            print("DATAFRAME_END")


        





    
    

    

##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
