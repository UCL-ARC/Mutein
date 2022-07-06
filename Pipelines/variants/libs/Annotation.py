


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
import CodingGene


class Annotation:
    def __init__(self,annotation_path,fsta):
        self.fasta = fsta
        self.annotation_path =annotation_path
        print("Loading genome...")
        #self.gencode = pd.read_table(self.annotation_path, comment="#", sep = "\t")        
        self.gencode = pd.read_table(self.annotation_path, comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])        
        #print("COLUMNS=", self.gencode.columns)
    
    def getCdsRegions(self,genes,gene_key,preferred_source=""):
        cds_dic = {} #per chme_gene per position/tuple NUC/AA
        gene_cds = {}
        print(self.gencode.columns)        
        have_gene = False
        cg_id = ""       
        coding_genes = []         
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
            cg = CodingGene.CodingGene(cg)
            chunks = []            
            print(cg)            
            for st,tpl in cds.items():                
                chme,strt,stp,strd,frm = tpl                                
                seq_chunk = self.fasta.getSeq(chme,int(strt),int(stp),strd=="-")
                print("seq chunk", strd=="-",seq_chunk)
                print("3 either side", self.fasta.getSeq(chme,int(strt)-3,int(stp)+3,strd=="-"))
                if strd=="-":
                    chunks.insert(0,[seq_chunk,tpl])
                else:
                    chunks.append([seq_chunk,tpl])
            print("---------------------------------------------")            
            seqaas = ""
            seqnucs = ""
            last_carry = ""
            for seq_chunk, tpl in chunks:                
                chme,strt,stp,strd,frm = tpl
                # frame includes push from last carry
                # frame is 3- for reverse
                if strd == "-" and int(frm) > 0:
                    frm = 3-int(frm)
                    frm = frm - len(last_carry)
                print("Chunk ends=",seq_chunk[:12],"...",seq_chunk[-12:])
                print("Carry over=",last_carry, "Phase=",frm)                                
                print(seq_chunk)
                seqnuc,seqaa,end_carry = codons.getAA(seq_chunk,last_carry,frm)
                print("Seq start=",seqnuc[0:3])
                print("Seq left=",divmod(len(seqnuc),3),end_carry)
                print(seqaa)                
                print(seqnuc)      
                print("try all phases")
                print(codons.getAA(seq_chunk,last_carry,0)[1])
                print(codons.getAA(seq_chunk,last_carry,1)[1])
                print(codons.getAA(seq_chunk,last_carry,2)[1])
                seqnucs += seqnuc
                seqaas += seqaa
                last_carry = end_carry
                cg.addChunk(int(strt),int(stp),strd=="+",seqnuc,seqaa)
                                

            seqnuc,seqaa,end_carry = codons.getAA(seqnuc)
            print("** NUC SEQ**")
            print(seqnuc)
            print("** AA SEQ**")
            print(seqaa)
            print("** NUC SEQS**")
            print(seqnucs)
            print("** AA SEQS**")
            print(seqaas)
            coding_genes.append(cg)
    
        return coding_genes

            

            

        
        
        #print("ids",len(dic_chm["ids"]),dic_chm["ids"])
        #df_vcf = pd.DataFrame.from_dict(dic_chm)
        #print(df_vcf)

