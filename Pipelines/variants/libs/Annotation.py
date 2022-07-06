


"""
RSA 30/6/22
We need to access the reference genome annotations  so need this data the gff3:
https://www.gencodegenes.org/human/
unzip: gzip -d gencode.v40.annotation.gff3.gz

Explanations: 
https://medium.com/intothegenomics/annotate-genes-and-genomic-coordinates-using-python-9259efa6ffc2
https://learn.gencore.bio.nyu.edu/ngs-file-formats/gff3-format/

"""

import pandas as pd
import codons


class Annotation:
    def __init__(self,annotation_path,fsta):
        self.fasta = fsta
        self.annotation_path =annotation_path
        print("Loading genome...")
        #self.gencode = pd.read_table(self.annotation_path, comment="#", sep = "\t")        
        self.gencode = pd.read_table(self.annotation_path, comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])        
        #print("COLUMNS=", self.gencode.columns)
    
    def getCdsRegions(self,genes,gene_key,preferred_source=""):
        gene_cds = {}
        print(self.gencode.columns)        
        have_gene = False
        cg_id = ""                
        for idx in self.gencode.index:
            chme = self.gencode["seqname"][idx]
            src = self.gencode["source"][idx]
            feat = self.gencode["feature"][idx]
            att = self.gencode["attribute"][idx]
            strt = self.gencode["start"][idx]
            stp = self.gencode["end"][idx]
            strd = self.gencode["strand"][idx]
            frm = self.gencode["frame"][idx]            
            atts = att.split(";")
            pref = True
            if preferred_source != "":            
                pref= preferred_source.upper()==src.upper()
                            
            if pref:
                if feat.upper() == "GENE":
                    have_gene=False
                    for at in atts:
                        if gene_key in at:
                            ats = at.split("=")
                            this_gene = ats[1].upper()
                            if this_gene in genes:
                                have_gene = True
                                cg_id = f"{chme}:{this_gene}".upper()                            
                                print(cg_id)
                                if cg_id not in gene_cds:
                                    gene_cds[cg_id] = {}
                                            
                if feat.upper() == "CDS" and have_gene:
                    # don't add it if the start date is the same as the last start date or one that already exists
                    # need to know if it is +ve or -ve
                    # and the phase
                    if strt not in gene_cds[cg_id]:
                        gene_cds[cg_id][strt] = [chme,strt,stp,strd,frm]
                        print(cg_id,src,strt,stp,strd,frm)                
                                        
        print("##################################################")
        
        for cg,cds in gene_cds.items():
            print(cg)
            seq = ""
            last_carry = ""
            for st,tpl in cds.items():
                print("Carry over=",last_carry)
                chme,strt,stp,strd,frm = tpl                            
                #def getSeq(self,chme,start,end,reverse=False):
                seq_chunk = self.fasta.getSeq(chme,strt,stp,strd=="-")
                seq += seq_chunk
                seqnuc,seqaa,end_carry = codons.getAA(seq_chunk,last_carry,frm)
                print("Seq start=",seqnuc[0:3])
                print("Seq left=",divmod(len(seq),3),end_carry)
                print(seqaa)
                last_carry = end_carry
                print("try all phases")
                print(codons.getAA(seq_chunk,last_carry,0))
                print(codons.getAA(seq_chunk,last_carry,1))
                print(codons.getAA(seq_chunk,last_carry,2))

            seqnuc,seqaa,end_carry = codons.getAA(seq)
            print("** NUC SEQ**")
            print(seqnuc)
            print("** AA SEQ**")
            print(seqaa)

            

        
        
        #print("ids",len(dic_chm["ids"]),dic_chm["ids"])
        #df_vcf = pd.DataFrame.from_dict(dic_chm)
        #print(df_vcf)

