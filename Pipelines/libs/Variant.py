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

    def includedInRange(self, starts, ends, offsets, chains):
        for s in range(len(starts)):
            start = int(starts[s])
            end = int(ends[s])
            offset = int(offsets[s])
            chain = chains[s]
            if self.residue >= int(start) and self.residue <= int(end):
                return True, offset, chain
        return False, 0, ""
