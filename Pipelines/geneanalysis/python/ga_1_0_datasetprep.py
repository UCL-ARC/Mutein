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
    dataset_path = Paths.Paths(data_dir, install_dir + "Pipelines/geneanalysis", dataset=dataset)
    
    import Pipelines.geneanalysis.python.ga_1_genestogene as ppla
    ppla.run_pipeline(args)

    # and we want only 1 batch for the stitching
    script_file = "libs/pipeline_qsubber.py"
    yaml_file = "geneanalysis/config/batch_gene_tasks.yml"
    bm = BatchMaker.BatchMaker(script_file, yaml_file)

    # load the list of the genes
    genes_fd = FileDf.FileDf(dataset_path.dataset_inputs + "genes_list.csv")    
    genes_csv = genes_fd.openDataFrame()
    genes = []
    print(genes_csv)
    for g in range(len(genes_csv.index)):
        gn = genes_csv["gene"][g]
        argsgn = args
        arglist = args[1]
        arglist += "@gene=" + gn
        argsgn[1] = arglist
        
        import Pipelines.geneanalysis.python.ga_2_genetoproteins as pplb
        print("Extracting pdbs for", gn)
        pplb.run_pipeline(argsgn)

        import Pipelines.geneanalysis.python.ga_2_genebackparams as pplc
        print("Extracting pdbs for", gn)
        exists = pplc.run_pipeline(argsgn)

        if exists:
            genes.append(gn)
            import Pipelines.geneanalysis.python.ga_2_genevarparams as ppld
            print("Extracting pdbs for", gn)
            ppld.run_pipeline(argsgn)

    genes_csv = FileDf.FileDic(dataset_path.dataset_inputs + "genes_pdb_list.csv",{})    
    for gn in genes:
        genes_csv.add("dataset",dataset)
        genes_csv.add("gene",gn)
        bm.addBatch(dataset, gn)

    genes_csv.saveAsDf()
            
    bm.printBatchScript(dataset_path.pipeline_path + "/foldx_"+ dataset + "_tasks.sh")

     
                            
    print("### COMPLETED dataset preparation ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
