
"""
RSA 30/6/22

class to test RefGenome
Find the gene seq here
https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606+AND+gene_exact:OXTR
"""

import pandas as pd
import Annotation
import FastaGenome
import FastaCds

genes = []
genes.append("ANG")
#genes.append("APOE")
#genes.append("APP")
#genes.append("NOTCH1")
#genes.append("WT1")
#this is a small -ve strand gene to test, not in our list
#genes.append("OXTR")

# Human genome 37
# Assemblies: https://www.ncbi.nlm.nih.gov/assembly/?term=hg19

# Human genome 38
# Assemblies https://www.ncbi.nlm.nih.gov/assembly/?term=hg38
fasta_path = "/home/rachel/UCL/github/MuteinData/data_genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
gff3_path = "/home/rachel/UCL/github/MuteinData/data_genome/Homo_sapiens.GRCh38.106.gff3"
pep_path = "/home/rachel/UCL/github/MuteinData/data_genome/Homo_sapiens.GRCh38.pep.all.fa"
cds_path = "/home/rachel/UCL/github/MuteinData/data_genome/Homo_sapiens.GRCh38.cds.all.fa"



fstCds = FastaCds.FastaCds(cds_path,pep_path)
seq = fstCds.getSeqPep("ANG")
print("***  TEST SEQ  ***")
print(seq)

fstGen = FastaGenome.FastaGenome(fasta_path)
anno = Annotation.Annotation(gff3_path,fstGen)


if True:
    
    cgs = anno.getCdsRegions(genes,"Name","ensembl_havana")
    print(cgs)
if False:
    anno.findMatchingRegion("WT1","MGSDVRDLNALLPAVPSLGGG","Name","ensembl_havana")
 



        