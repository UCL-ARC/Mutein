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
is_dir=False
# 37
pep_path = "/home/rachel/UCL/github/MuteinData/data_genome/37/Homo_sapiens.GRCh37.pep.all.fa"
cds_path = "/home/rachel/UCL/github/MuteinData/data_genome/37/Homo_sapiens.GRCh37.cds.all.fa"
fasta_path = "/home/rachel/UCL/github/MuteinData/data_genome/37/assembly/"
gff3_path = "/home/rachel/UCL/github/MuteinData/data_genome/37/Homo_sapiens.GRCh37.87.gff3"
is_dir=True



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
#genes = ["APOE"]
#genes = ["ADAM29"]
#genes = ["ANG"]
#genes = ["APP"]
#genes = ["NOTCH1"]
#genes = ["WT1"] #no matches
#genes = ["OXTR"]

# test cases
genes = ["CR2"] #forward strand
#genes = ["AURKA"] #small and negative


matches = []
not_matches = []

# Find all the genes that we are able to deal with in the first version
for gene in genes:
    print(gene)
    seq = cm_search.seqFromUniProt(gene)
    found,gene_name,seq,tr = fstCds.getSeqDetailsPep(gene.upper(),seq)
    if found:
        matches.append([gene,tr,seq])
        #print("Gene found",gene,f"\n{seq}")
    else:
        not_matches.append(gene)

    
print("------------------")
print("Matches")
print(matches)
print("------------------")
if len(not_matches) > 0:
    print("!!! Not Matches !!!")
    print(not_matches)
    print("------------------")

fstGen = FastaGenome.FastaGenome(fasta_path,is_dir)
anno = Annotation.Annotation(gff3_path,fstGen)

print("########  Creating the CDS #############")
for gene,tr,seq in matches:
    coding_gene = anno.getCdsRegions(tr)
    print(coding_gene.getString(loglevel=0))
    #print(coding_gene.aas)
    if not seq==coding_gene.aas:
        #print("Correct=",seq)
        print(gene,"No match, Found=",f"\n{coding_gene.aas}")
        print(gene,"Orig=",f"\n{seq}")
    else:
        print(gene,"FOUND A MATCH :-)")
        if gene.upper() == "CR2":  #positive
            """
            207640071	A	G	1 207640071 A G	0.000279263	CR2	CCDS31007.1	r.448a>g	c.259A>G	p.K87E
            207642037	G	A	1 207642037 G A	0.000383403	CR2	CCDS31007.1	r.800g>a	c.611G>A	p.S204N
            207646227	G	A	1 207646227 G A	0.000379219	CR2	CCDS31007.1	r.1870g>a	c.1681G>A	p.E561K
            207647013	C	T	1 207647013 C T	0.000374065	CR2	CCDS31007.1	r.2291c>u	c.2102C>T	p.S701F
            """
            retsA= coding_gene.getVariant(207640071,"G",True)            
            retsB= coding_gene.getVariant(207642037,"A",True)            
            retsC= coding_gene.getVariant(207646227,"A",True)
            retsD= coding_gene.getVariant(207647013,"T",True)
            print(retsA)
            print(retsB)
            print(retsC)
            print(retsD)


    
    



