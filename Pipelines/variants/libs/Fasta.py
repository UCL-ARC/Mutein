"""
RSA 5/7/22
------------------------
Class to manage gene and chromosome position

"""
from Bio import SeqIO

class Fasta:
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.seqs = {}
        

    def getSeq(self,chme,start,end,reverse=False):
        cid = f"CHR{chme}"        
        if cid not in self.seqs:
            for seq_record in SeqIO.parse(self.fasta_file, "fasta"):
                if seq_record.id.upper() == cid:
                    #return seq_record.seq[start-1:end]
                    self.seqs[cid] = seq_record.seq
                    break
                    
        if cid in self.seqs:            
            if reverse:
                rev_seq = self.flipSeq(self.seqs[cid][start-1:end][::-1])
                return rev_seq
            else:
                seq = self.seqs[cid][start-1:end]
                return seq
        else:
            return ""

    def flipSeq(self,seq):    
        rev_seq = ""
        for n in seq:            
            fn = self.flipNucleotide(n)
            rev_seq += fn
        return rev_seq#[::-1]

    def flipNucleotide(self,nuc):    
        if nuc.upper() == "C":
            return "G"
        elif nuc.upper() == "G":
            return "C"
        elif nuc.upper() == "A":
            return "T"
        elif nuc.upper() == "T":
            return "A"
        else:
            return "?"

    
                
