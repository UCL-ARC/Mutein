"""
RSA 5/7/22
------------------------
Class to manage gene and chromosome position

"""
from Bio import SeqIO

class FastaGenome:
    def __init__(self, fasta_file, is_dir=False):
        print("Creating genome file data from fasta")
        self.chmes = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
        self.seqs = {}
        if is_dir:
            for chm in self.chmes:
                cid = f"CHR{chm}"       
                chm_file = f"{fasta_file}/chr{chm}.fna" 
                for seq_record in SeqIO.parse(chm_file, "fasta"):
                    self.seqs[cid] = seq_record.seq
                    #print(str(seq_record.seq))
        else:               
            for seq_record in SeqIO.parse(fasta_file, "fasta"):
                cid = str(seq_record.id).upper()
                if cid not in self.seqs:                                                      
                    self.seqs[cid] = str(seq_record.seq)
                    #print("Storing",str(cid))
                    
                
    def getSeq(self,chme,start,end,reverse=False):
        cid = f"CHR{chme}"                                    
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

    
                
