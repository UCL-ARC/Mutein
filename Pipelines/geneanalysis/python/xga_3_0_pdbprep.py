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
    gene = argus.arg("gene")
    pdb = argus.arg("pdb")
    dataset_path = Paths.Paths(data_dir, install_dir + "Pipelines/geneanalysis", dataset=dataset,gene=gene,pdb=pdb)
                    
    import Pipelines.geneanalysis.python.mod_pdbparams as pplc
    print("Back params for", pdb)
    exists = pplc.run_pipeline(args)
    
    import Pipelines.geneanalysis.python.mod_pdbvparams as ppld
    print("Variant params for", pdb)
    ppld.run_pipeline(args)
                                     
    print("### COMPLETED pdb preparation ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys
    globals()["run_pipeline"](sys.argv)
