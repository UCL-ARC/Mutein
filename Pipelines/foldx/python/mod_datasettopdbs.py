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
import FileDf


def run_pipeline(args):
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    dataset_path = Paths.Paths(
        data_dir, install_dir, dataset=dataset
    )

    # We want batches for prep and tasks
    #script_file = "libs/pipeline_qsubber.py"
    #pdbs_yaml_file = "geneanalysis/config/batch_gene_1_pdbs.yml"
    #pdbs_bm = BatchMaker.BatchMaker(script_file, pdbs_yaml_file)
    #rep_yaml_file = "geneanalysis/config/batch_gene_2_rep.yml"
    #rep_bm = BatchMaker.BatchMaker(script_file, rep_yaml_file)
    #prep_yaml_file = "geneanalysis/config/batch_gene_3_prep.yml"
    #prep_bm = BatchMaker.BatchMaker(script_file, prep_yaml_file)
    #tasks_yaml_file = "geneanalysis/config/batch_gene_4_tasks.yml"
    #tasks_bm = BatchMaker.BatchMaker(script_file, tasks_yaml_file)

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

        import mod_genetoproteins as pplb

        print("Extracting pdbs for", gn)
        pdbs = pplb.run_pipeline(argsgn)
        if pdbs > 0:
            genes.append(gn)
        """
        import Pipelines.geneanalysis.python.ga_2_genebackparams as pplc
        print("Extracting pdbs for", gn)
        exists = pplc.run_pipeline(argsgn)

        if exists:
            genes.append(gn)
            import Pipelines.geneanalysis.python.ga_2_genevarparams as ppld
            print("Extracting pdbs for", gn)
            ppld.run_pipeline(argsgn)

        """
    genes_csv = FileDf.FileDic(dataset_path.dataset_inputs + "genes_pdb_list.csv", {})
    for gn in genes:
        genes_csv.add("dataset", dataset)
        genes_csv.add("gene", gn)
        #tasks_bm.addBatch(dataset, gn)
        #prep_bm.addBatch(dataset, gn)
        #pdbs_bm.addBatch(dataset, gn)
        #rep_bm.addBatch(dataset, gn)

    genes_csv.saveAsDf()

    """
    pdbs_bm.printBatchScript(
        dataset_path.pipeline_path + "/foldx_" + dataset + "_1_pdbs.sh"
    )
    rep_bm.printBatchScript(
        dataset_path.pipeline_path + "/foldx_" + dataset + "_2_rep.sh"
    )
    prep_bm.printBatchScript(
        dataset_path.pipeline_path + "/foldx_" + dataset + "_3_prep.sh"
    )
    tasks_bm.printBatchScript(
        dataset_path.pipeline_path + "/foldx_" + dataset + "_4_tasks.sh"
    )
    """

    print("### COMPLETED dataset preparation ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
