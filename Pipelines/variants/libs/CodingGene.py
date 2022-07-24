"""
RSA 6/7/22

This class is to try t make it easier to work with CDS regions of genes
In particular, to find variants from chromosome positions
given that the variants could be on the forward or reverse strands
"""

import codons


class CodingGene:
    def __init__(self, transcript):
        self.transcript = transcript    
        self.chunks = {}
        self.aas = ""
        self.nucs = ""
                                
    def getString(self,loglevel=0):
        ret = self.transcript + "\nPositions"
        for pos,tpl in self.chunks.items():
            ret += "\n" + str(tpl[0]) + ":" + str(tpl[1])
        if loglevel>1:
            ret += "\nNucSeq:"+self.nucs
        if loglevel>0:            
            ret += "\nAASeq:"+self.aas
        return ret

    def addChunk(self,start,end,fwd,nucleotides,aminos):
        # we add a position with start and end
        #if fwd:
        self.chunks[int(start)] = [int(start),int(end),nucleotides,aminos,fwd]
        #else:
            # it is kind of horrible that it needs to be reversed back
            # it is the only way to get the position correctly
            #self.chunks[int(start)] = [int(start),int(end),nucleotides[::-1],aminos[::-1]]            
    
    def getAminoAcid(self,chmepos):
        start,end,nucs,aas = self.getTuple(chmepos)
        relative_pos = chmepos-start
        div,rem = divmod(len(relative_pos),3)
        return aas[div+1]

    def getNucleotide(self,chmepos):
        start,end,nucs,aas,fwd = self.getTuple(chmepos)
        relative_pos = chmepos-start
        return nucs[relative_pos]

    def getVariant(self,chmepos,mut,fwd=True):
        start,end,nucs,aas,fwd = self.getTuple(chmepos)
        print("tup ret",start,end,nucs,aas,fwd)
        relative_pos = chmepos-start
        orig_nuc = nucs[relative_pos]
        print("orig nuc",orig_nuc)
        div,rem = divmod(relative_pos,3)
        print("div rem",div,rem)
        orig_aa = aas[div+1]
        if rem == 0:
            triple = nucs[relative_pos:relative_pos+3]
        elif rem == 1:
            triple = nucs[relative_pos-1:relative_pos+2]
        elif rem == 2:
            triple = nucs[relative_pos-3:relative_pos+1]
        if not fwd:
            triple = triple[::-1]
        #def getAA(triple_seq,append_first="", phase=0):
        print("triple",triple)
        new_aa = codons.getAA(triple)
        print("new aa",new_aa)
        #nuc, nuc_mut, aa, aa_mu, aa_num
        return orig_aa,new_aa
        
    def getTuple(self,chmepos):
        for pos,tpl in self.chunks.items(): #tpl=[int(start),int(end),nucleotides,aminos,fwd]
            start,end,nucleotides,aminos,fwd = tpl
            print("Seek",chmepos,"from",start,end,nucleotides,aminos,fwd)
            if int(chmepos) >= start and int(chmepos) <= end:
                return start,end,nucleotides,aminos,fwd
        return [0,0,"","",True]

    



