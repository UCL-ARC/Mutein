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
    dataset_path = Paths.Paths(
        data_dir, install_dir + "Pipelines/geneanalysis", dataset=dataset
    )
    genes_list = []
    """
    "gene_name"	"n_syn"	"n_mis"	"n_non"	"n_spl"	"n_ind"	"wmis_cv"	"wnon_cv"	"wspl_cv"	"wind_cv"	"pmis_cv"	"ptrunc_cv"	"pallsubs_cv"	"pind_cv"	"qmis_cv"	"qtrunc_cv"	"qallsubs_cv"	"pglobal_cv"	"qglobal_cv"
    "NOTCH1"	404	2700	772	360	1208	3.98846230912502	21.7266421919647	21.7266421919647	17.6547328391658	0	0	0	8.07429581976078e-68	0	0	0	0	0
    """
    # 2.) First get list of genes:
    if len(genes_list) == 0:
        genes_file = dataset_path.dataset_inputs + "/genes.txt"
        with open(genes_file, "r") as fr:
            lines = fr.readlines()
            # lines = ["","FAT1"]
            for l in range(1, len(lines)):
                # for l in range(1,2):
                line = lines[l]
                cols = line.split("\t")
                gene = cols[0].upper()
                if gene[0] == '"':
                    gene = gene[1:]
                if gene[-1] == '"':
                    gene = gene[:-1]
                genes_list.append(gene)

    # 1.) First prepare the variants file
    genes_variant_file = dataset_path.dataset_inputs + "genes_variants.txt"
    gene_variant_dic = genestovariants.extractVariantsFromFile(genes_variant_file)
    genes = []
    for gene in genes_list:
        gene_path = Paths.Paths(
            data_dir, install_dir + "Pipelines/geneanalysis", dataset=dataset, gene=gene
        )
        accession = genetoprotein.accession_from_bioservices(gene.upper())
        if len(accession) > 1:
            seq = genetoprotein.sequence_from_bioservices(accession)
            seq_lines = seq.split("\n")
            wholeseq = ""
            for s in range(1, len(seq_lines)):
                sl = str(seq_lines[s].strip())
                wholeseq += sl
            gn = Gene.Gene(gene, accession, wholeseq)
            genes.append(gn)  # main repository for data we are creating in function
            script_file = "libs/pipeline_qsubber.py"
            yaml_file = "geneanalysis/config/batch_pdb.yml"
            bm = BatchMaker.BatchMaker(script_file, yaml_file)
            # CREATE the variants for the gene
            vrs = gene_variant_dic[gene.upper()]
            for i in range(len(vrs["bases"])):
                bs = vrs["bases"][i]
                p1 = vrs["prev_aa"][i]
                rs = vrs["residue"][i]
                na = vrs["new_aa"][i]
                vr = vrs["variant"][i]
                # print(bs,p1,rs,na,vr)
                vrnt = Variant.Variant(gene, vr, bs)
                gn.addVariant(vrnt)
            dfv = gn.getVariantsDataFrame()
            dfv.to_csv(gene_path.gene_inputs + "/variants.csv", index=False)

    # make a list of the genes
    genes_csv = FileDf.FileDic(dataset_path.dataset_inputs + "genes_list.csv", {})
    for gn in genes_list:
        genes_csv.add("dataset", dataset)
        genes_csv.add("gene", gn)
    genes_csv.saveAsDf()

    print("### COMPLETED genes to gene ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
