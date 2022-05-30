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

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)
import Paths
import Arguments
import BatchMaker
import Gene
import Variant
import genetoprotein
import genestovariants
import SwissModel
import UniProt
import PdbRunner
import FileDf


def run_pipeline(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")
    sys.path.append(install_dir)
    sys.path.append(install_dir + "/Pipelines")
    sys.path.append(install_dir + "/Pipelines/libs")
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene","")
    pdb = argus.arg("pdb","")
    runs = argus.arg("runs")
        
    runner = None
    if pdb != "":
        import Pipelines.geneanalysis.python.ga_pdb_runner as runner        
    elif gene != "":
        import Pipelines.geneanalysis.python.ga_gene_runner as runner
    else:
        import Pipelines.geneanalysis.python.ga_dataset_runner as runner

    if "a" in runs:
        import Pipelines.geneanalysis.python.ga_dataset_runner as runnerd
        print("Mutein: Preparing genes")
        runnerd.prepareGenes(args)
    if "b" in runs:
        print("Mutein: Preparing pdbs")
        runner.preparePdbs(args)
    if "c" in runs:
        print("Mutein: Repairing pdbs")
        runner.repairPdbs(args)
    if "d" in runs:
        print("Mutein: Making background param file")
        runner.makeParams(args)
    if "e" in runs:
        print("Mutein: Making variant param file")
        runner.makeVparams(args)
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
        import Pipelines.geneanalysis.python.ga_gene_runner as runnerj
        print("Mutein: Aggregating gene tasks")
        runnerj.aggGene(args)
                                                                 
    print("### COMPLETED Mutein script ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
