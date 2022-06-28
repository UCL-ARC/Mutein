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


def run_pipeline(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset", "")
    gene = argus.arg("gene")
    pdb = argus.arg("pdb")    
    for pdb in [pdb]:
        pdb_path = Paths.Paths(
                    data_dir,
                    install_dir,
                    dataset=dataset,
                    gene=gene,
                    pdb=pdb,
                )
        
        pdburl = genetoprotein.getPDBLink(pdb)
        biopdb = genetoprotein.retrievePdbStructure(pdburl, pdb, pdb_path.pdb_inputs + pdb + ".pdb")
        
    print("### COMPLETED gene to proteins pipeline ###")
    print("MUTEIN SCRIPT ENDED")
    return [pdb]


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
