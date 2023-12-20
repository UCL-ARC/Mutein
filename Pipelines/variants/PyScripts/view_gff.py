
import pandas as pd
genes = []
genes.append("ANG")
genes.append("APOE")
genes.append("APP")
genes.append("NOTCH1")

vcf_path = "/home/rachel/UCL/github/MuteinData/data_genome/gencode.v40.annotation.gff3"
fasta_path = "/home/rachel/UCL/github/MuteinData/data_genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
print("Loading genome...")
gencode = pd.read_table(vcf_path, comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])        
print("...loaded genome")
print(gencode.columns)
gene_cds = {}
for gene in genes:
    print("Searching for",gene)
    for idx in gencode.index:
        chme = gencode["seqname"][idx]
        feat = gencode["feature"][idx]
        att = gencode["attribute"][idx]
        strt = gencode["start"][idx]
        stp = gencode["end"][idx]
        atts = att.split(";")
        cg_id = (chme+":").upper()
        if feat == "CDS":
            for at in atts:
                if "gene_name" in at:
                    ats = at.split("=")
                    if ats[1].upper() == gene.upper():                                                     
                        cg_id += gene.upper()
                        print(cg_id)
                        if cg_id not in gene_cds:
                            gene_cds[cg_id] = []
                        gene_cds[cg_id].append([strt,stp])
                        
print(gene_cds)

                            
dic_chm = {}
dic_chm["ids"] = []
dic_chm["cds"] = []
dic_chm["length"] = []
dic_chm["seq"] = []

print("##################################################")
for cg,cds in gene_cds.items():
    print(cg,cds)
    dic_chm["ids"].append(cg)
    ss = ""
    ln = 0
    for cd in cds:
        print(cd)
        ss += f"{cd[0]}:{cd[1]}|"
        ln += cd[1]-cd[0] + 1
    dic_chm["cds"].append(ss)
    dic_chm["length"].append(ln)
    dic_chm["seq"].append("")


print("ids",len(dic_chm["ids"]),dic_chm["ids"])
df_vcf = pd.DataFrame.from_dict(dic_chm)
print(df_vcf)
