
import pandas as pd
from Bio import SeqIO

fasta_path = "/home/rachel/UCL/github/MuteinData/data_genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"


for seq_record in SeqIO.parse(fasta_path, "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

