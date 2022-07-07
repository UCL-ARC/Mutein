"""
RSA 7/7/22
------------------------
Class to manage amino acid protein sequence from genome

"""
from Bio import SeqIO

class FastaCds:
    def __init__(self, fasta_file_cds, fasta_file_pep):
        self.fasta_file_cds = fasta_file_cds
        self.fasta_file_pep = fasta_file_pep
        self.seqs = {}
        
    def getSeq(self,gene_name):
        seqs = []
        gid = f"gene_symbol:{gene_name.upper() }"                
        for seq_record in SeqIO.parse(self.fasta_file_cds, "fasta"):
            if gid in seq_record.description:                
                seqs.append([seq_record.seq,seq_record.description])
        return seqs

    def getSeqPep(self,gene_name):
        seqs = []
        gid = f"gene_symbol:{gene_name.upper()} description"                
        for seq_record in SeqIO.parse(self.fasta_file_pep, "fasta"):
            if gid in seq_record.description:                
                seqs.append([seq_record.seq,seq_record.description])
        return seqs
                
                    
                
