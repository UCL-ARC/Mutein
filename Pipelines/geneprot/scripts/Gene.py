"""
RSA 12/4/22
------------------------
Class to manage the data associated with a gene

"""
import os
import pandas as pd

class Gene:
    def __init__(self, gene,seq):
        self.gene = gene
        self.sequence = seq
        self.variants = {}
        self.pdbs = {}
    
    def addVariant(self,variant):
        if variant.variant not in self.variants:
            self.variants[variant.variant] = variant                    
    
    def addPdb(self,pdb):
        if pdb.pdb not in self.pdbs:
            self.pdbs[pdb.pdb] = pdb

    def getMatchingVarantPdb(self):
        v_p = []
        for v in self.variants:
            for p in self.pdbs:
                if p.matchesResidue(v.residue):
                    v_p.append([v,p])
        return v_p
    
    def getVariantCandidatesDataFrame(self):
        dic_variants = {}
        dic_variants["gene"] = []
        dic_variants["variant"] = []
        dic_variants["residue"] = []
        dic_variants["bases"] = []
        dic_variants["candidates"] = []        
        
        for vrcod,vr in self.variants.items():
            candidate = ""
            print(vr.variant)
            for pdbcod,pdb in self.pdbs.items():
                if pdb.matchesResidue(vr.residue):
                    candidate += pdbcod + " "
            dic_variants["residue"].append(vr.residue)
            dic_variants["gene"].append(self.gene)
            dic_variants["candidates"].append(candidate)
            dic_variants["variant"].append(vr.variant)
            dic_variants["bases"].append(vr.bases)
        
        pdbs_df = pd.DataFrame.from_dict(dic_variants)
        pdbs_df = pdbs_df.sort_values(by='residue', ascending=True)
        return pdbs_df
        
    