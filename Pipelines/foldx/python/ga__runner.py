"""
RSA 11.4.22

Script to take a file in the format given to me by Michael Hall and then produce a list of proteins

We use the bioservices python library to access databases
https://pypi.org/project/bioservices/

"""
import os
import sys
import yaml
import pandas as pd

import _helper
import Arguments

def run_pipeline(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset","")
    gene = argus.arg("gene", "")
    pdb = argus.arg("pdb", "")
    runs = argus.arg("runs")

    runner = None
    if pdb != "":
        import ga_pdb_runner as runner
    elif gene != "":
        import ga_gene_runner as runner
    else:
        import ga_dataset_runner as runner
    
    import mod_cleaning as cleaner

    if "a" in runs:        
        print("Mutein: Preparing genes")
        runner.prepareGenes(args)
    if "b" in runs:
        print("Mutein: Preparing pdbs")
        runner.preparePdbs(args)
    if "c" in runs:
        print("Mutein: Repairing pdbs")
        success = runner.repairPdbs(args)
        if success:
            cleaner.cleanPdbThruputs(args)
    if "d" in runs:
        print("Mutein: Making background param file")
        runner.makeParams(args)
        cleaner.cleanPdbThruputs(args)
    if "e" in runs:
        print("Mutein: Making variant param file")
        runner.makeVparams(args)
        cleaner.cleanPdbThruputs(args)
    if "f" in runs:
        print("Mutein: Running background tasks")
        runner.runTasks(args)        
    if "g" in runs:
        print("Mutein: Running variant tasks")
        runner.runVtasks(args)        
    if "h" in runs:
        print("Mutein: Aggregating background tasks")
        runner.aggTasks(args)
    if "i" in runs:
        print("Mutein: Aggregating variant tasks")
        runner.aggVtasks(args)
    if "j" in runs:        
        print("Mutein: Aggregating gene tasks")
        runner.aggGene(args)
    if "k" in runs:        
        print("Mutein: Aggregating dataset tasks")
        runner.aggGenes(args)
    if "x" in runs:        
        print("Mutein: Cleaning data")
        cleaner.cleanGene(args)

    print("### COMPLETED Mutein script ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys
    globals()["run_pipeline"](sys.argv)