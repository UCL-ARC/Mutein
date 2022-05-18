"""
RSA 11/5/22
------------------------
Class to manage coverage of a pdb to a gene

"""
import os
from urllib.request import urlretrieve
import Bio.PDB as bio
import warnings

warnings.filterwarnings("ignore")  # sadly because of the annoying bioython warnings


class PdbCoverage:
    def __init__(self,struc,seq):
        self.struc = struc
        self.seq = seq
    
    def getCoverage(self,minfrag=-1):
        minfrag = -1 #no matter what is passed in, it only works with whole fragments at the moment
        # PPBuilder is C-N and CAPPBuilder is CA-CA
        ppb = bio.CaPPBuilder()                        
        has_match = False
        peptides = []
        try:
            peptides = ppb.build_peptides(self.struc)            
        except:
            return []
        
        # we will return a tuple that has chain, pdb start, pdb end, gene start, gene end
        segments = []

        for pp in peptides:
            seq_one = str(pp.get_sequence())            
            resis = pp.get_ca_list()[0]
            chain = resis.parent.get_parent().id
            # https://biopython.org/docs/1.75/api/Bio.PDB.Atom.html
            # get_full_id(self): Return the full id of the atom.
            # The full id of an atom is the tuple (structure id, model id, chain id, residue id, atom name, altloc).
            residue_num = resis.get_full_id()[3][1]            
            residue_end = residue_num + len(seq_one)-1            

            
            if minfrag == -1 or minfrag > len(seq_one):
                fragment = len(seq_one)
            else:
                fragment = minfrag
            
            start = 0            
            while start < len(seq_one)-fragment+1:
                a_match = False
                end = start+fragment            
                seq_frag = seq_one[start:end]            
                gene_start = self.seq.find(seq_frag)+1
                if gene_start > 0:
                    a_match=True
                    end = start+fragment+1
                    matches=True
                    while end < len(seq_one)+1 and matches:
                        seq_frag = seq_one[start:end]            
                        gene_start_new = self.seq.find(seq_frag)+1
                        if gene_start_new > 0:
                            matches=True
                            end += 1
                            gene_start=gene_start_new
                        else:
                            matches=False
                            end -=1   
                            seq_frag = seq_frag[:-1]                                                         
                if a_match:
                    resis = pp.get_ca_list()[start]
                    residue_num = pp.get_ca_list()[start].get_full_id()[3][1]            
                    chain = resis.parent.get_parent().id
                    #residue_num += start
                    residue_end = residue_num + len(seq_frag)-1
                    gene_end = gene_start + len(seq_frag)-1
                    coverage = round(len(seq_frag)/len(self.seq),4)
                    segment = [chain,residue_num,residue_end,gene_start,gene_end,coverage]
                    print(segment)
                    segments.append(segment)
                start=end

        return segments
            