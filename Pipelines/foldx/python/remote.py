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
import sys

import _helper
import Arguments
import Paths
import BatchStatus
from os import listdir
from os.path import isfile, join
import os


def checkResult(onefile):
    if exists(onefile):
        timeA = pathlib.Path(onefile).stat().st_mtime
        return True, datetime.fromtimestamp(timeA).strftime("%d%b%y-%H:%M")
    return False, ""


def checkResults(ddg, bm):
    if exists(ddg):
        timeA = pathlib.Path(ddg).stat().st_mtime
        print(
            "DDG Background file was created at",
            datetime.fromtimestamp(timeA).strftime("%d%b%y-%H:%M"),
        )
    else:
        print("!!!DDG Background file does not exist")
    if exists(bm):
        timeB = pathlib.Path(bm).stat().st_mtime
        print(
            "Variant buildmodel file was created at",
            datetime.fromtimestamp(timeB).strftime("%d%b%y-%H:%M"),
        )
    else:
        print("!!!Variant buildmodel file does not exist")
    
def run_pipeline(args):
    now = datetime.now()
    current_time = now.strftime("%d%b%y-%H.%H.%S")
    print("Mutein remote script:", current_time)

    ##############################################
    mode = args[1]
    pattern = args[2]
    WorkDir = args[3]
    DataDir = args[4]
    InstallDir = args[5]
    PipelineDir = args[6]

    dataset_gene_pdb = pattern.split(":")
    dataset, gene, pdb = (
            dataset_gene_pdb[0],
            dataset_gene_pdb[1],
            dataset_gene_pdb[2],
        )
    
    # if it is a comma delim gene list then split it more
    genes = gene.split(",")
    for gene in genes:        
        print("Mode=", mode)
        print("Pattern=", pattern)
        print(gene)
        if mode == "GENES":
            batch_stat = BatchStatus.BatchStatus(DataDir, PipelineDir,dataset,"","")
            batch_stat.createReport()       
            break 
        elif mode == "PDBS":
            batch_stat = BatchStatus.BatchStatus(DataDir, PipelineDir,dataset,gene,"")
            batch_stat.createReport()              
        elif mode == "PDB":
            batch_stat = BatchStatus.BatchStatus(DataDir, PipelineDir,dataset,gene,pdb)
            batch_stat.createReport()              
        elif mode == "DSINCOMPLETE":
            batch_stat = BatchStatus.BatchStatus(DataDir, PipelineDir,dataset,"","")
            genelist = batch_stat.getGenes()    
            for gene in genelist:
                pdblist = batch_stat.getGenePdbs(gene)    
                for pdb in pdblist:          
                    batch_stat.createUntasksForPdb(gene,pdb)        
            break
        elif mode == "DSINCOMPLETEPDB":
            batch_stat = BatchStatus.BatchStatus(DataDir, PipelineDir,dataset,"","")
            genelist = batch_stat.getGenes()    
            for gene in genelist:
                batch_stat.makeMissingPdbs(gene)                
            break
        elif mode == "GENEINCOMPLETEPDB":        
            batch_stat = BatchStatus.BatchStatus(DataDir, PipelineDir,dataset,gene,"")
            batch_stat.makeMissingPdbs(gene)        
        elif mode == "PDBINCOMPLETE":
            dataset_gene_pdb = pattern.split(":")
            dataset, gene, pdb = (
                dataset_gene_pdb[0],
                dataset_gene_pdb[1],
                dataset_gene_pdb[2],
            )
            batch_stat = BatchStatus.BatchStatus(DataDir, PipelineDir,dataset,gene,pdb)
            batch_stat.createUntasksForPdb(gene,pdb)                
            break
        elif mode == "GENEINCOMPLETE":            
            batch_stat = BatchStatus.BatchStatus(DataDir, PipelineDir,dataset,gene,"")
            pdblist = batch_stat.getGenePdbs(gene)#
            #We also need to create the genes level incomplete file
            path = Paths.Paths(DataDir, PipelineDir, dataset=dataset, gene=gene)
            filename_incomplete_b = path.inputs + "params_background_incomplete.txt"
            filename_incomplete_v = path.inputs + "params_variants_incomplete.txt"
            count = 0                
            with open(filename_incomplete_b, "w") as fw_back:
                with open(filename_incomplete_v, "w") as fw_var:
                    for pdb in pdblist:          
                        batch_stat.createUntasksForPdb(gene,pdb,True,fw_back,fw_var,count)    
                        count += 1
        elif mode == "PDB_BACK":            
            dataset_gene_pdb = pattern.split(":")
            dataset, gene, pdb = (
                dataset_gene_pdb[0],
                dataset_gene_pdb[1],
                dataset_gene_pdb[2],
            )            
            path = Paths.Paths(DataDir, PipelineDir, dataset=dataset, gene=gene,pdb=pdb)
            if "GENE" in mode or pdb == "":            
                filename = path.outputs + "ddg_variants.csv"
            else:
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
        elif mode == "GENE_BACK" or mode=="GENE_VAR":            
            dataset_gene_pdb = pattern.split(":")
            dataset, gene, pdb = (
                dataset_gene_pdb[0],
                dataset_gene_pdb[1],
                dataset_gene_pdb[2],
            )
            ds_path = Paths.Paths(DataDir,PipelineDir,dataset=dataset)
            ge_path = Paths.Paths(DataDir,PipelineDir,dataset=dataset,gene=gene)
            if "VAR" in mode:                                
                filename = ds_path.outputs + f"{gene.lower()}_ddg_variants.csv"
            else:                
                filename = ds_path.outputs + f"{gene.lower()}_ddg_background.csv"            
            mexists, time = checkResult(filename)
            if not mexists:
                if "VAR" in mode:                
                    filename = ge_path.outputs + "ddg_variants.csv"
                else:
                    filename = ge_path.outputs + "ddg_background.csv"
            if mexists:
                print(filename,"contains results")
                print("DATAFRAME_START")
                with open(filename, "r") as fr:
                    lines = fr.readlines()
                    for line in lines:
                        print(line.strip())
                print("DATAFRAME_END")
            else:
                print(filename "does not exist")

        elif mode == "PDB_BM":
            dataset_gene_pdb = pattern.split(":")
            dataset, gene, pdb = (
                dataset_gene_pdb[0],
                dataset_gene_pdb[1],
                dataset_gene_pdb[2],
            )
            path = Paths.Paths(DataDir, PipelineDir, dataset=dataset, gene=gene, pdb=pdb)
            filename = path.outputs + "ddg_variants"                        
            mexists, time = checkResult(filename)
            if mexists:
                print(filename)
                print("DATAFRAME_START")
                with open(filename, "r") as fr:
                    lines = fr.readlines()
                    for line in lines:
                        print(line.strip())
                print("DATAFRAME_END")        
        elif mode == "COVERAGE":
            dataset_gene_pdb = pattern.split(":")
            dataset, gene, pdb = (
                dataset_gene_pdb[0],
                dataset_gene_pdb[1],
                dataset_gene_pdb[2],
            )
            path = Paths.Paths(DataDir, PipelineDir, dataset=dataset, gene=gene)
            filename = path.outputs + "Coverage_all.csv"
            mexists, time = checkResult(filename)
            if mexists:
                print("DATAFRAME_START")
                with open(filename, "r") as fr:
                    lines = fr.readlines()
                    for line in lines:
                        print(line.strip())
                print("DATAFRAME_END")
        elif mode == "DS_GENES":
            dataset_gene_pdb = pattern.split(":")
            dataset, gene, pdb = (
                dataset_gene_pdb[0],
                dataset_gene_pdb[1],
                dataset_gene_pdb[2],
            )
            path = Paths.Paths(DataDir, PipelineDir, dataset=dataset)
            filename = path.inputs + "genes_pdb_list.csv"
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
