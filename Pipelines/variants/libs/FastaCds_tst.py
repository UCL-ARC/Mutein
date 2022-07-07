"""
RSA 7/7/22
------------------------
Class to test FastaCds regions

"""
import FastaCds

# 38
pep_path = "/home/rachel/UCL/github/MuteinData/data_genome/Homo_sapiens.GRCh38.pep.all.fa"
cds_path = "/home/rachel/UCL/github/MuteinData/data_genome/Homo_sapiens.GRCh38.cds.all.fa"

# 37
pep_path = "/home/rachel/UCL/github/MuteinData/data_genome/37/Homo_sapiens.GRCh37.pep.all.fa"
cds_path = "/home/rachel/UCL/github/MuteinData/data_genome/37/Homo_sapiens.GRCh37.cds.all.fa"

fstCds = FastaCds.FastaCds(cds_path,pep_path)

"""
seqs = fstCds.getSeqPep("ANG")
print("***  ANG  ***")
for seq in seqs:
    print("------")
    print(len(seq),seq)

seqs = fstCds.getSeqPep("APOE")
print("***  APOE  ***")
for seq in seqs:
    print("------")
    print(len(seq),seq)

"""
seqs = fstCds.getSeqPep("WT1")
print("***  WT1  ***")
for seq in seqs:
    print("------")
    print(len(seq[0]))
    print(seq[1])
    print(seq[0])