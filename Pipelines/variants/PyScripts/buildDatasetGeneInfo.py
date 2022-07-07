"""
RSA 6/7/21

This file takes a gene list as input (from the results of Rob's pipeline)
This file is simple a list of genes

I then build it up to show the gene sequence, length, and coding regions so that I can identify variants

Headers

Chme,Gene,Accession,UniprotSeq,UniprotLength,CodingRegions,GeneNucSeq,GeneNucLength,GeneAASeq,GeneAALength

"""
import os
import sys
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]
lib_path = "/".join(dirs) + "/libs"
sys.path.append(lib_path)
import Chme
import Fasta
import Annotation

# I may want to build this file up incrementally
#I start with opening whatever I have got and I I am re-running a particular bit I write over it

hardcoded_gene_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/gene_list.txt"
genes = []

# Get the list of genes we are working with
with open(hardcoded_gene_path, mode='r') as org:
    lines = org.readlines()
    for line in lines:
        gene = line.strip()
        if len(gene)>1:
            genes.append(gene)

hardcoded_ouput_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/gene_annotated.csv"
##)Open up the file if it exists and get the data in it currently - it is in order and MUST have genes

##)
fasta_path = "/home/rachel/UCL/github/MuteinData/data_genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
fst = Fasta.Fasta(fasta_path)

anno = Annotation.Annotation("/home/rachel/UCL/github/MuteinData/data_genome/Homo_sapiens.GRCh38.106.gff3",fst)
cgs = anno.getCdsRegions(genes,"Name","ensembl_havana")
print(cgs)