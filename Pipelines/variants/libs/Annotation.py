


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
        self.gencode = pd.read_table(self.annotation_path, comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'],low_memory=False)        
        #print("COLUMNS=", self.gencode.columns)
    
    
    def getCdsRegions(self,transcript,loglevel=0):        
        cdg = CodingGene.CodingGene(transcript)
        trans_starts = {}
        seqaas = ""
        seqnucs = ""
        #print(self.gencode.columns)                
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
                                                                                                                    
            if feat.upper() == "CDS":
                # don't add it if the start date is the same as the last start date or one that already exists
                # need to know if it is +ve or -ve
                # and the phase
                #print(atts)
                for at in atts:
                    if "transcript:" in at:
                        ats = at.split(":")
                        this_trn = ats[1].upper()
                        if this_trn == transcript:                                                        
                            trans_starts[strt] = [chme,strt,stp,strd,frm]
                            print(this_trn,chme,strt,stp,strd,frm)
                            
                                                            
        if loglevel > 1:
            print("##################################################")        
        
        chunks = []   
        for strt,tpl in trans_starts.items():
            print(strt,tpl)                                                           
            chme,strt,stp,strd,frm = tpl
            if loglevel > 1:
                print("Seek chme", chme,strt,stp,strd=="-")                       
            seq_chunk = self.fasta.getSeq(chme,int(strt),int(stp),strd=="-")
            if loglevel > 1:
                print("seq chunk", strd=="-",seq_chunk)
                print("3 either side", self.fasta.getSeq(chme,int(strt)-3,int(stp)+3,strd=="-"))
            if strd=="-":
                chunks.insert(0,[seq_chunk,tpl])
            else:
                chunks.append([seq_chunk,tpl])
        
        if loglevel > 1:
            print("---------------------------------------------")                        
        last_carry = ""
        for seq_chunk, tpl in chunks:                
            if loglevel > -1:
                print("## Chunk",tpl,"last_carry=",last_carry)
            chme,strt,stp,strd,frm = tpl
            # frame includes push from last carry
            # frame is 3- for reverse            
            if False:
                if strd == "-" and int(frm) > 0:
                    frm = 3-int(frm)
                    frm = frm - len(last_carry)
            # match off the carry and the frameshift.
            #If forward and a carry of 2 then there is no frameshift
            if int(frm) == 2 and len(last_carry) == 1 and strd == "+":
                frm = 0
            if loglevel > -1:
                print("Chunk ends=",seq_chunk[:12],"...",seq_chunk[-5:])                                
            seqnuc,seqaa,end_carry = codons.getAA(seq_chunk,last_carry,int(frm))
            if loglevel > 1:
                print("Seq ends=",seqnuc[:12],"...",seqnuc[-12:])                                                
                print("-------------Sequences---------------")
                print(seqaa)     
                print("----------------------------")
                print(seqnuc)      
                print("-------------Try all phases---------------")
                print(codons.getAA(seq_chunk,last_carry,0)[1])
                print(codons.getAA(seq_chunk,last_carry,1)[1])
                print(codons.getAA(seq_chunk,last_carry,2)[1])
            seqnucs += seqnuc                
            seqaas += seqaa
            last_carry = end_carry                
            cdg.addChunk(int(strt),int(stp),strd=="+",seqnuc,seqaa)                                    
        
        # the whole thing ######################
        seqnuc,seqaa,end_carry = codons.getAA(seqnuc)
        if loglevel > 1:
            print("** NUC SEQ**")
            print(seqnuc)
            print("** AA SEQ**")
            print(seqaa)
            print("** NUC SEQS**")
            print(seqnucs)
            print("** AA SEQS**")
            print(seqaas)
        cdg.nucs = seqnucs
        if seqaas[-1]=="X":
            cdg.aas = seqaas[:-1]
        else:
            cdg.aas = seqaas
                
        return cdg
        
    def findMatchingRegion(self,gene,seq_in,gene_key,preferred_source=""):                        
        have_gene = False        
        gene_cds = {}              
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
                            if this_gene == gene:
                                have_gene = True
                                cg_id = f"{chme}:{this_gene}".upper()                            
                                print(cg_id)
                                if cg_id not in gene_cds:
                                    gene_cds[cg_id] = {}
                                            
                if feat.upper() == "CDS" and have_gene:                    
                    if strt not in gene_cds[cg_id]:
                        gene_cds[cg_id][strt] = [chme,strt,stp,strd,frm]
                        print(cg_id,src,strt,stp,strd,frm)                
                                        
        print("##################################################")        
        for cg,cds in gene_cds.items():            
            chunks = []            
            print(cg)            
            for st,tpl in cds.items():                
                chme,strt,stp,strd,frm = tpl                                
                seq_chunk = self.fasta.getSeq(chme,int(strt),int(stp),strd=="-")                                
                if strd=="-":
                    chunks.insert(0,[seq_chunk,tpl])
                else:
                    chunks.append([seq_chunk,tpl])
            print("---------------------------------------------")            
            seqaas = ""
            seqnucs = ""
            last_carry = ""
            for seq_chunk, tpl in chunks:                
                print("#################  New chunk",tpl,"last_carry=",last_carry)
                chme,strt,stp,strd,frm = tpl                                
                print("Chunk ends=",seq_chunk[:12],"...",seq_chunk[-12:])                                                
                print("-------------Try all phases---------------")
                ph1 = codons.getAA(seq_chunk,last_carry,0)[1]
                ph2 = codons.getAA(seq_chunk,last_carry,1)[1]
                ph3 = codons.getAA(seq_chunk,last_carry,2)[1]
                if ph1[2:] in seq_in:
                    print("*******************  Found phase 1 in seq!", ph1)
                elif ph2[2:] in seq_in:
                    print("*******************  Found phase 2 in seq!", ph2)
                elif ph3[2:] in seq_in:
                    print("*******************  Found phase 3 in seq!", ph3)
                elif seq_in[2:] in ph1:
                    print("*******************  Found seq in phase 1!", ph1)
                elif seq_in[2:] in ph2:
                    print("*******************  Found seq in phase 2!", ph2)
                elif seq_in[2:] in ph3:
                    print("*******************  Found seq in phase 3!", ph3)
                
                                

            
            
    
        

            

            

    
    
    #print("ids",len(dic_chm["ids"]),dic_chm["ids"])
    #df_vcf = pd.DataFrame.from_dict(dic_chm)
    #print(df_vcf)

