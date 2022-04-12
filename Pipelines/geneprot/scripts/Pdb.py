"""
RSA 12/4/22
------------------------
Class to manage the data associated with a pdb file specifically for a gene

"""
import os

class Pdb:
    def __init__(self,gene,pdb,chain,segment_start,segment_end,method, resolution):
        self.gene = gene
        self.pdb = pdb
        self.chain = chain
        self.segment_start = int(segment_start)
        self.segment_end = int(segment_end)
        self.method = method
        self.resolution = resolution

    def matchesResidue(self,residue):
        if int(residue) >= self.segment_start and int(residue) <= self.segment_end:
            return True
        else:
            return False
