"""
RSA 6/7/21

This file takes a gene list as input (from the results of Rob's pipeline)
This file is simple a list of genes

I then build it up to show the gene sequence, length, and coding regions so that I can identify variants

Headers

Chme,Gene,Accession,UniprotSeq,UniprotLength,CodingRegions,GeneNucSeq,GeneNucLength,GeneAASeq,GeneAALength

"""
# I may want to build this file up incrementally
#I start with opening whatever I have got and I I am re-running a particular bit I write over it

hardcoded_gene_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/gene_list"
genes = []

# Get the list of genes we are working with
with open(hardcoded_gene_path, mode='r') as org:
    lines = org.readlines()
    gene = lines.strip()
    if len(gene)>1:
        genes.append()

hardcoded_ouput_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/gene_annotated.csv"
#Open up the file if it exists and get the data in it currently - it is in order and MUST have genes


