"""
RSA 30/6/22
Use the packakge scikit-allel
(pip install scikit-allel)
https://scikit-allel.readthedocs.io/en/stable/
Some instructions here:
http://alimanfoo.github.io/2017/06/14/read-vcf.html




"""

import sys
import os


dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]
lib_path = "/".join(dirs) + "/libs"
sys.path.append(lib_path)
import Chme


import allel
import pandas as pd

## INPUTS ##
hardcoded_vcf_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/keogh2018.vcf"
hardcoded_gene_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/gene_list"
hardcoded_interval_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/intervals.list"
hardcoded_organism_id = "/home/rachel/UCL/github/MuteinData/vcf_keogh/organism_id.txt"
## OUPUTS ##
hardcoded_out_path = "/home/rachel/UCL/github/MuteinData/dataset_keogh/genes_variants_inputs.csv"
print("################ CONSTRUCT GENE DICTIONARY ###############")
chrms = []
starts = []
ends = []
genes = []
organism_id = ""
with open(hardcoded_organism_id, mode='r') as org:
    lines = org.readlines()
    organism_id = lines[0].strip()
with open(hardcoded_gene_path, mode='r') as gr:
    lines = gr.readlines()
    for lin in lines:
        line = lin.strip()
        genes.append(line)
with open(hardcoded_interval_path, mode='r') as gr:
    lines = gr.readlines()
    for lin in lines:
        line = lin.strip()        
        c1 = line.split(":")
        chrms.append(c1[0])
        c2 = c1[1].split("-")
        starts.append(c2[0])
        ends.append(c2[1])
chromos = {}
if not (len(chrms) == len(genes) and len(starts) == len(ends) and len(ends) == len(chrms)):
    raise Exception("Problem with input data lengths")

for i in range(len(chrms)):
    ch = chrms[i]
    gn = genes[i]
    st = starts[i]
    en = ends[i]
    if ch not in chromos:
        chm = Chme.Chme(ch,organism_id)
        chromos[ch] = chm
    chm = chromos[ch]
    chm.addGene(gn,int(st),int(en))        

chm = chromos["chr1"]
chm.findGene(123456)         

print("################ PARSE VCF FILE ###############")
with open(hardcoded_vcf_path, mode='r') as org:
    print(org.readlines()[:50])
dic_vcf = {}
callset = allel.read_vcf(hardcoded_vcf_path,fields="variants/*")
for keyx,arr in callset.items():    
    if "variants" in keyx:
        key = keyx[9:]
        print(key)        
        shp = arr.shape
        if len(shp) > 1:
            x,y = shp[0],shp[1]            
            vals = []            
            for i in range(x):
                val = ""
                for j in range(y):                    
                    val += str(arr[i,j])+":"
                vals.append(val[:-1])
            dic_vcf[key] = vals        
        else:            
            dic_vcf[key] = arr

# make it a dataframe
df_vcf = pd.DataFrame.from_dict(dic_vcf)
# now we need to add the gene column on to it and some other manipulations
#mouse,sample,chr,pos,from_base,to_base,gene,impact,variant
#MD5690,MD5690f,2,26462159,G,A,Notch1,nonsense,p.Q2047*
"""
'impact == "missense"'
variant,from_base,to_base
"""
genes = []
impacts = []
tos = []
for idx in df_vcf.index:
    pos = df_vcf["POS"][idx]
    chm = df_vcf["CHROM"][idx]
    mss = df_vcf["is_snp"][idx]
    to = df_vcf["ALT"][idx]
    chme = chromos[chm]
    gn = chme.findGene(int(pos))    
    genes.append(gn)        
    if len(gn) > 1:
        if mss:
            impacts.append("missense")
            tos.append(to[0])
            # now map the nucleotides to the residues
            seq = chme.seqFromUniProt(gn)                    
        else:
            impacts.append("other")
            tos.append(to)
    else:
        impacts.append("no_gene")
        tos.append(to)
df_vcf["gene"] = genes    
df_vcf["impact"] = impacts
df_vcf["from_base"] = df_vcf["REF"]
df_vcf["to_base"] = tos

    
print("Save VCF dataframe to",hardcoded_out_path)
df_vcf.to_csv(hardcoded_out_path,index=False)
