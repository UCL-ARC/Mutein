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
    gene_path = Paths.Paths(data_dir, install_dir + "Pipelines/geneanalysis", dataset=dataset,gene=gene)
            
            
    import Pipelines.geneanalysis.python.mod_genetoproteins as pplb
    print("Extracting pdbs for", gene)
    pplb.run_pipeline(args)

    import Pipelines.geneanalysis.python.mod_geneparams as pplc
    print("Extracting pdbs for", gene)
    exists = pplc.run_pipeline(args)

    if exists:        
        import Pipelines.geneanalysis.python.mod_genevparams as ppld
        print("Extracting pdbs for", gene)
        ppld.run_pipeline(args)
                                            
    print("### COMPLETED gene preparation ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys
    globals()["run_pipeline"](sys.argv)
