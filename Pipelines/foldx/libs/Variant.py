"""
RSA 12/4/22
------------------------
Class to manage the data associated with a gene's variant

"""
import os


class Variant:
    def __init__(self, gene, variant, bases):
        self.gene = gene
        from_aa = variant[:1]
        to_aa = variant[-1:]
        residue = variant[1:-1]
        self.residue = int(residue)
        self.from_aa = from_aa
        self.to_aa = to_aa
        self.variant = variant
        self.bases = bases

    def includedInRange(self, segments):
        for chain, residue_num, residue_end, gene_start, gene_end, coverage in segments:
            if self.residue >= int(gene_start) and self.residue <= int(gene_end):
                offset = gene_start - residue_num
                res_start = self.residue - offset
                return True, chain, res_start
        return False, 0, ""
