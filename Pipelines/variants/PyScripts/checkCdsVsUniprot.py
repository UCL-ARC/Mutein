"""
RSA 7/7/21

Given a list of genes, I want to find the uniprot sequence and then all the CDS transcirptions
If any of them match, hooray!
Otherwise flag it

"""
import os
import sys
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]
lib_path = "/".join(dirs) + "/libs"
sys.path.append(lib_path)
import UniProt
import FastaCds
import FastaGenome
import Annotation


genes = []
#genes.append("ANG")
#genes.append("APOE")
#genes.append("APP")
#genes.append("NOTCH1")
#genes.append("WT1")
#this is a small -ve strand gene to test, not in our list
#genes.append("OXTR")
pep_path = "/home/rachel/UCL/github/MuteinData/data_genome/38/106/Homo_sapiens.GRCh38.pep.all.fa"
cds_path = "/home/rachel/UCL/github/MuteinData/data_genome/38/106/Homo_sapiens.GRCh38.cds.all.fa"
fasta_path = "/home/rachel/UCL/github/MuteinData/data_genome/38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
gff3_path = "/home/rachel/UCL/github/MuteinData/data_genome/38/106/Homo_sapiens.GRCh38.106.gff3"

fstCds = FastaCds.FastaCds(cds_path,pep_path)
cm_search = UniProt.UniProt(9606)

hardcoded_gene_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/gene_list.txt"
# Get the list of genes we are working with
with open(hardcoded_gene_path, mode='r') as org:
    lines = org.readlines()
    for line in lines:
        gene = line.strip()
        if len(gene)>1:
            genes.append(gene)


##### Or replace with chosen genes ###
genes = ["NOTCH1"]

matches = []
not_matches = []

# Find all the genes that we are able to deal with in the first version
for gene in genes:
    print(gene)
    seq = cm_search.seqFromUniProt(gene)
    found = False
    seqs = fstCds.getSeqPep(gene.upper())
    for seqi in seqs:
        seqo = seqi[0]
        if seqo == seq: 
            if not found:
                found = True    
                ids = seqi[1].split(" ")
                trs = ""
                for ide in ids:
                    if "transcript:" in ide:
                        trss = ide.split(":")
                        trs = trss[1].split(".")
                        tr = trs[0]
                print(trs)
                matches.append([gene,seq,tr])
    
    if not found:                    
        not_matches.append(gene)

    
print("------------------")
print("Matches")
print(matches)
print("------------------")
if len(not_matches) > 0:
    print("!!! Not Matches !!!")
    print(not_matches)
    print("------------------")

fstGen = FastaGenome.FastaGenome(fasta_path)
anno = Annotation.Annotation(gff3_path,fstGen)

print("########  Creating the CDS #############")
transcripts = []
for gene, seq,transcript in matches:
    print(gene, transcript)
    transcripts.append(transcript)
    #print(fstGen.getSeq(1,500,520,False))
    #print(fstGen.getSeq(19,44906625,44906667,False))
    #print(fstGen.getSeq(19,44907760,44907952,False))
    #print(fstGen.getSeq(19,44908533,44909250,False))
for transcript in transcripts:
    #maybe I want to get the regions out of anno and only after that ask for fasta anf the AFTER that get the protein seqeunce
    anno.getCdsRegions(transcripts)
    



