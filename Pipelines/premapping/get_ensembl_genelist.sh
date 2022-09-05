#!/bin/bash

# manually download the list of human gene positions and names from ensembl
# using biomart as explained in the readme.md in this folder
# save it as mart_export.txt.gz
#the file should contain 7 columns exactly matching:
#geneid,geneid version,chromosome,gene start, gene end,gene name,gene synonym

#uncompress and rename
gunzip ensembl_gene_list/mart_export.txt.gz
mv ensembl_gene_list/mart_export.txt ensembl_gene_list/all_genes.tsv

#make a version of the list lacking synonyms
cat ensembl_gene_list/all_genes.tsv | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' \
| sort -u > all_genes_no_synonyms.tsv

#extract gene positions for keogh2018 where gene_list contains the 102 gene names
#from supplementary table 1
awk 'NR == FNR {genes[$0]} NR > FNR && ($6 in genes)' \
datasets/keogh2018/gene_list ensembl_gene_list/all_genes_no_synonyms.tsv  \
| grep -v CHR_ | sort -u > datasets/keogh2018/keogh2018_gene_positions.tsv

#print only the position and full chromosome name
awk '{print "chr"$3":"$4"-"$5}' datasets/keogh2018/keogh2018_gene_positions.tsv \
> datasets/keogh2018/intervals.list