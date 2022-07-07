"""
RSA 6/7/22

This class is to try t make it easier to work with CDS regions of genes
In particular, to find variants from chromosome positions
given that the variants could be on the forward or reverse strands
"""

import codons


class CodingGene:
    def __init__(self, gene_name):
        self.gene = gene_name        
        self.chunks = {}
        self.aas = ""
        self.nucs = ""
                                
    def getString(self):
        ret = self.gene + "/n"
        for pos,tpl in self.chunks.items():
            ret += str(tpl)
        return ret



    def addChunk(self,start,end,fwd,nucleotides,aminos):
        # we add a position with start and end
        if fwd:
            self.chunks[int(start)] = [int(start),int(end),nucleotides,aminos]
        else:
            # it is kind of horrible that it needs to be reversed back
            # it is the only way to get the position correctly
            self.chunks[int(start)] = [int(start),int(end),nucleotides[::-1],aminos[::-1]]            
    
    def getAminoAcid(self,chmepos):
        start,end,nucs,aas = self.getTuple(chmepos)
        relative_pos = chmepos-start
        div,rem = divmod(len(relative_pos),3)
        return aas[div+1]

    def getNucleotide(self,chmepos):
        start,end,nucs,aas = self.getTuple(chmepos)
        relative_pos = chmepos-start
        return nucs[relative_pos]

    def getVariant(self,chmepos,mut,fwd=True):
        start,end,nucs,aas = self.getTuple(chmepos)
        relative_pos = chmepos-start
        div,rem = divmod(len(relative_pos),3)
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
        new_aa = codons.getAA(triple)
        return orig_aa,new_aa
        
    def getTuple(self,chmepos):
        for pos,tpl in self.chunks.items():
            if int(chmepos) >= pos:
                return tpl
        return [0,0,"",""]

    



