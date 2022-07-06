
"""
RSA 30/6/22

class to test RefGenome
"""

import pandas as pd
import Annotation
import Fasta

genes = []
#genes.append("ANG")
genes.append("APOE")
#genes.append("APP")
#genes.append("NOTCH1")


fasta_path = "/home/rachel/UCL/github/MuteinData/data_genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
fst = Fasta.Fasta(fasta_path)

anno = Annotation.Annotation("/home/rachel/UCL/github/MuteinData/data_genome/Homo_sapiens.GRCh38.106.gff3",fst)
anno.getCdsRegions(genes,"Name","ensembl_havana")

#anno = Annotation.Annotation("/home/rachel/UCL/github/MuteinData/data_genome/gencode.v40.annotation.gff3")
#anno.getCdsRegions(genes,"gene_name")

        