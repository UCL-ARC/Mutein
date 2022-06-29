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
    dataset = argus.arg("dataset","")
    gene = argus.arg("gene","")
    genes_list = []
    if gene != "":
        genes_list = [gene]
            
    dataset_path = Paths.Paths(data_dir, install_dir, dataset=dataset,readonly=False)
    genes_variant_file = dataset_path.dataset_inputs + "genes_variants_inputs.csv"    
    
    """
    "gene_name"	"n_syn"	"n_mis"	"n_non"	"n_spl"	"n_ind"	"wmis_cv"	"wnon_cv"	"wspl_cv"	"wind_cv"	"pmis_cv"	"ptrunc_cv"	"pallsubs_cv"	"pind_cv"	"qmis_cv"	"qtrunc_cv"	"qallsubs_cv"	"pglobal_cv"	"qglobal_cv"
    "NOTCH1"	404	2700	772	360	1208	3.98846230912502	21.7266421919647	21.7266421919647	17.6547328391658	0	0	0	8.07429581976078e-68	0	0	0	0	0
    """
    # The file contains the variants and the ist of genes
    
    genes_df = FileDf.FileDf(genes_variant_file, sep=",", header=True).openDataFrame()

    # this contains the organism id from uniport, eg human=9606 and mouse=10090
    organism_id_path = dataset_path.dataset_inputs + "organism_id.txt"
    with open(organism_id_path, "r") as fr:
        lines = fr.readlines()
        organism_id = lines[0].strip()
    
    # 2.) First get list of genes:
    acc_no = {}
    acc_gene = {}
    if len(genes_list) == 0:                
        genes_list = genes_df["gene"].unique()

    # 1.) First prepare the variants file    
    gene_variant_dic = genestovariants.extractVariantsFromFile(genes_df)
    genes = []
    for gene in genes_list:
        gene_path = Paths.Paths(
            data_dir, install_dir, dataset=dataset, gene=gene,,readonly=False
        )
        accessions = genetoprotein.accession_from_bioservices(gene.upper(),organism_id,True)
        if len(accessions) == 0:
            accessions = genetoprotein.accession_from_bioservices(gene.upper(),organism_id,False)
        if len(accessions) > 0:
            for accession in accessions:
                seq = genetoprotein.sequence_from_bioservices(accession)
                seq_lines = seq.split("\n")
                wholeseq = ""
                for s in range(1, len(seq_lines)):
                    sl = str(seq_lines[s].strip())
                    wholeseq += sl
                gn = Gene.Gene(gene, accession, wholeseq)
                genes.append(gn)  # main repository for data we are creating in function                
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
                acc_gene[gn.accession] = gene
                acc_no[gn.accession] = gn.numVariants()
                dfv = gn.getVariantsDataFrame()
                dfv.to_csv(gene_path.gene_inputs + "/variants.csv", index=False)
                
    # make a list of the genes
    genes_csv = FileDf.FileDic(dataset_path.dataset_inputs + "genes_list.csv", {})
    genes_csv_all = FileDf.FileDic(dataset_path.dataset_inputs + "genes_list_all.csv", {})
    
    used_genes = []
    
    for w in sorted(acc_no, key=acc_no.get, reverse=True):
        acc,no = w, acc_no[w]
        gn = acc_gene[acc]
        genes_csv_all.add("dataset", dataset)
        genes_csv_all.add("gene", gn)
        genes_csv_all.add("organism", organism_id)
        genes_csv_all.add("accession", acc)
        genes_csv_all.add("no_vars", no)

        if gn not in used_genes:
            genes_csv.add("dataset", dataset)
            genes_csv.add("gene", gn)
            genes_csv.add("organism", organism_id)
            genes_csv.add("accession", acc)
            genes_csv.add("no_vars", no)
            used_genes.append(gn)
        
            

    #for gn,no in genes_dic.items():
    #    genes_csv.add("dataset", dataset)
    #    genes_csv.add("gene", gn)
    #    genes_csv.add("no_vars", no)
    genes_csv.saveAsDf()
    genes_df = genes_csv_all.saveAsDf()
    for gene in genes_list:
        genes_one = genes_df.query('gene == "'+gene+'"')        
        gene_path = Paths.Paths(data_dir, install_dir, dataset=dataset, gene=gene,readonly=False)
        one_path = gene_path.gene_inputs + "accessions.csv"
        genes_one.to_csv(one_path, index=False)
        

    #print("### COMPLETED genes to gene ###")
    #print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)
