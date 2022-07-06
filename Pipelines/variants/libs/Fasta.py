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
            seq = self.seqs[cid][start-1:end]
            if reverse:
                return seq[::-1]
            else:
                return seq
        else:
            return ""

    
                
