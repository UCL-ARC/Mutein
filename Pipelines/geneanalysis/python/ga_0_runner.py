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
        runner = PdbRunner.PdbRunner()
    elif gene != "":
        runner = GeneRunner.GeneRunner()
    else:
        runner = DatasetRunner.DatasetRunner()

    if "a" in runs:
        runner.prepareGenes()
    if "b" in runs:
        runner.preparePdbs()
    if "c" in runs:
        runner.repairPdbs()
    if "d" in runs:
        runner.makeParams()
    if "e" in runs:
        runner.makeVparams()
    if "f" in runs:
        runner.runTasks()
    if "g" in runs:
        runner.runVtasks()
    if "h" in runs:
        runner.aggTasks()
    if "i" in runs:
        runner.aggVtasks()
    if "j" in runs:
        runner.aggTasks()
    if "k" in runs:
        runner.aggGene()

    

    
    
    
            
    
                                 
    print("### COMPLETED dataset preparation ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
