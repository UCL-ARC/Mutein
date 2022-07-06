
"""
RSA 30/6/22

class to test RefGenome
Find the gene seq here
https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606+AND+gene_exact:OXTR
"""

import pandas as pd
import Annotation
import Fasta

genes = []
#genes.append("ANG")
#genes.append("APOE")
#genes.append("APP")
genes.append("NOTCH1")

#this is a small -ve strand gene to test, not in our list
#genes.append("OXTR")



fasta_path = "/home/rachel/UCL/github/MuteinData/data_genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
fst = Fasta.Fasta(fasta_path)

anno = Annotation.Annotation("/home/rachel/UCL/github/MuteinData/data_genome/Homo_sapiens.GRCh38.106.gff3",fst)
cgs = anno.getCdsRegions(genes,"Name","ensembl_havana")
print(cgs)

#anno = Annotation.Annotation("/home/rachel/UCL/github/MuteinData/data_genome/gencode.v40.annotation.gff3")
#anno.getCdsRegions(genes,"gene_name")

        