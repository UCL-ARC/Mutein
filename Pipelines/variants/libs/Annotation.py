


"""
RSA 30/6/22
We need to access the reference genome annotations  so need this data the gff3:
https://www.gencodegenes.org/human/
unzip: gzip -d gencode.v40.annotation.gff3.gz

Explanation: https://medium.com/intothegenomics/annotate-genes-and-genomic-coordinates-using-python-9259efa6ffc2
"""

import pandas as pd


class Annotation:
    def __init__(self,annotation_path):
        self.annotation_path =annotation_path
        print("Loading genome...")
        self.gencode = pd.read_table(self.annotation_path, comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])        
        print(self.gencode.columns)
    
    def getCdsRegions(self,genes):
        gene_cds = {}
        for gene in genes:
            print("Searching for",gene)
            for idx in self.gencode.index:
                chme = self.gencode["seqname"][idx]
                feat = self.gencode["feature"][idx]
                att = self.gencode["attribute"][idx]
                strt = self.gencode["start"][idx]
                stp = self.gencode["end"][idx]
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
