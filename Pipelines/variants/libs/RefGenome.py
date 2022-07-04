


"""
RSA 30/6/22
We need to access the reference genome annotations  so need this data the gff3:
https://www.gencodegenes.org/human/
unzip: gzip -d gencode.v40.annotation.gff3.gz

Explanation: https://medium.com/intothegenomics/annotate-genes-and-genomic-coordinates-using-python-9259efa6ffc2
"""

import pandas as pd


class RefGenome:
    def __init__(self):
        self.path = "/home/rachel/UCL/github/MuteinData/vcf/gencode.v40.annotation.gff3"
        gencode = pd.read_table(self.path, comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])        
        #gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)
        chms = []
        genes = []
        starts = []
        stops = []
        for idx in gencode.index:
            chme = gencode["seqname"][idx]
            feat = gencode["feature"][idx]
            att = gencode["attribute"][idx]
            strt = gencode["start"][idx]
            atts = att.split(";")            
            if "start_codon" in feat:
                for at in atts:
                    chms.append(chme)
                    starts.append(strt)
                    if "gene_id" in at:
                        genes.append(at)
            if "stop_codon" in feat:
                stops.append(strt)

        dic_chm = {}
        dic_chm["chme"] = chme
        dic_chm["gene"] = genes
        dic_chm["start"] = starts
        dic_chm["stop"] = stops

        print("CHME",len(chme))
        print("GENE",len(genes))
        print("START",len(starts))
        print("STOP",len(stops))

        df_vcf = pd.DataFrame.from_dict(dic_chm)
        print(df_vcf)
                                        

        