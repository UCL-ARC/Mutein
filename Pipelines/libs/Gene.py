"""
RSA 12/4/22
------------------------
Class to manage the data associated with a gene

"""
import os
import pandas as pd

import Pdb


class Gene:
    def __init__(self, gene, accession, seq):
        self.gene = gene
        self.accession = accession
        self.sequence = seq
        self.variants = {}
        self.pdbs = {}

    def addVariant(self, variant):
        if variant.variant not in self.variants:
            self.variants[variant.variant] = variant

    def getPdb(self, pdbcode, method, res):
        if pdbcode not in self.pdbs:
            pdb = Pdb.Pdb(self.gene, pdbcode, method, res)
            self.pdbs[pdbcode] = pdb
        else:
            pdb = self.pdbs[pdbcode]
        return pdb
    
    def addPdb(self, pdb):        
        self.pdbs[pdb.pdbcode] = pdb
        
    def getMatchingVarantPdb(self):
        v_p = []
        for v in self.variants:
            for p in self.pdbs:
                if p.matchesResidue(v.residue):
                    v_p.append([v, p])
        return v_p

    def getVariantCandidatesDataFrame(self):
        dic_variants = {}
        dic_variants["gene"] = []
        dic_variants["accession"] = []
        dic_variants["variant"] = []
        dic_variants["residue"] = []
        dic_variants["bases"] = []
        dic_variants["candidates"] = []

        for vrcod, vr in self.variants.items():
            candidate = ""
            # print(vr.variant)
            for pdbcod, pdb in self.pdbs.items():
                if pdb.matchesResidue(vr.residue):
                    candidate += pdbcod + " "
            dic_variants["residue"].append(vr.residue)
            dic_variants["gene"].append(self.gene)
            dic_variants["accession"].append(self.accession)
            dic_variants["candidates"].append(candidate)
            dic_variants["variant"].append(vr.variant)
            dic_variants["bases"].append(vr.bases)

        pdbs_df = pd.DataFrame.from_dict(dic_variants)
        pdbs_df = pdbs_df.sort_values(by="residue", ascending=True)
        return pdbs_df
    
    def getVariantsDataFrame(self):
        dic_variants = {}
        dic_variants["gene"] = []
        dic_variants["accession"] = []
        dic_variants["variant"] = []
        dic_variants["gene_no"] = []
        dic_variants["bases"] = []
        
        for vrcod, vr in self.variants.items():            
            dic_variants["gene"].append(self.gene)
            dic_variants["accession"].append(self.accession)            
            dic_variants["variant"].append(vr.variant)
            dic_variants["gene_no"].append(vr.residue)                                    
            dic_variants["bases"].append(vr.bases)        

        vars_df = pd.DataFrame.from_dict(dic_variants)
        vars_df = vars_df.sort_values(by="gene_no", ascending=True)
        return vars_df


    def getPdbCoverageDataFrame(self):
        dic_coverage = {}
        dic_coverage["gene"] = []
        dic_coverage["accession"] = []
        dic_coverage["pdb"] = []
        dic_coverage["method"] = []
        dic_coverage["resolution"] = []
        dic_coverage["coverage"] = []

        for pdbcod, pdb in self.pdbs.items():
            dic_coverage["gene"].append(self.gene)
            dic_coverage["accession"].append(self.accession)
            dic_coverage["pdb"].append(pdbcod)
            dic_coverage["method"].append(pdb.getMethod())
            dic_coverage["resolution"].append(pdb.resolution)
            dic_coverage["coverage"].append(pdb.getSegments())
        pdbs_df = pd.DataFrame.from_dict(dic_coverage)
        return pdbs_df

    def getSinglePdbCoverageDataFrame(self, pdb):
        dic_coverage = {}
        dic_coverage["gene"] = []
        dic_coverage["accession"] = []
        dic_coverage["pdb"] = []
        dic_coverage["chain"] = []
        dic_coverage["method"] = []
        dic_coverage["resolution"] = []
        dic_coverage["pdb_start"] = []
        dic_coverage["pdb_end"] = []
        dic_coverage["offset"] = []
        dic_coverage["gene_start"] = []
        dic_coverage["gene_end"] = []

        # dic_coverage["coverage"].append(pdb.getSegments())
        # segs += seg_start + ":" + seg_end + ":" + seg_off + ":" + self.chain + " "
        segments = pdb.getSegments().split(" ")
        for segment in segments:
            segs = segment.split(":")
            if len(segs) >= 4:
                gene_start = int(segs[0])
                gene_end = int(segs[1])
                seg_off = int(segs[2])
                seg_chain = segs[3]
                pdb_start = gene_start - seg_off
                pdb_end = gene_end - seg_off

                dic_coverage["gene"].append(self.gene)
                dic_coverage["accession"].append(self.accession)
                dic_coverage["pdb"].append(pdb.pdbcode)
                dic_coverage["chain"].append(seg_chain)
                dic_coverage["method"].append(pdb.getMethod())
                dic_coverage["resolution"].append(pdb.resolution)
                dic_coverage["pdb_start"].append(pdb_start)
                dic_coverage["pdb_end"].append(pdb_end)
                dic_coverage["offset"].append(seg_off)
                dic_coverage["gene_start"].append(gene_start)
                dic_coverage["gene_end"].append(gene_end)
        pdbs_df = pd.DataFrame.from_dict(dic_coverage)
        return pdbs_df

    def getPdbVariantCoverageDataFrame(self, pdb):
        dic_coverage = {}
        dic_coverage["gene"] = []
        dic_coverage["accession"] = []
        dic_coverage["pdb"] = []
        dic_coverage["variant"] = []
        dic_coverage["chain"] = []
        dic_coverage["gene_no"] = []
        dic_coverage["gene_mut"] = []        
        dic_coverage["pdb_residue"] = []
        dic_coverage["pdb_mut"] = []

        for vrcod, vr in self.variants.items():
            included, chain, pdb_residue = vr.includedInRange(pdb.segments)                            
            if included:
                dic_coverage["gene"].append(self.gene)
                dic_coverage["accession"].append(self.accession)
                dic_coverage["pdb"].append(pdb.pdbcode)
                dic_coverage["variant"].append("ANY")
                dic_coverage["chain"].append(chain)
                dic_coverage["gene_no"].append(vr.residue)
                dic_coverage["gene_mut"].append(vrcod)                
                dic_coverage["pdb_residue"].append(int(pdb_residue))
                dic_coverage["pdb_mut"].append(vrcod[0]+str(pdb_residue)+vrcod[-1]
                )

        pdbs_df = pd.DataFrame.from_dict(dic_coverage)
        pdbs_df = pdbs_df.sort_values(by="gene_no", ascending=True)
        return pdbs_df

    def getDatasetGenesPdbsDataFrame(self, dataset, genes):
        dic_coverage = {}
        dic_coverage["dataset"] = []
        dic_coverage["gene"] = []
        dic_coverage["accession"] = []
        dic_coverage["variants"] = []
        dic_coverage["length"] = []
        dic_coverage["pdbs"] = []

        for gene in genes:
            dic_coverage["dataset"].append(dataset)
            dic_coverage["gene"].append(gene.gene)
            dic_coverage["accession"].append(gene.accession)
            dic_coverage["variants"].append(len(gene.variants.items()))
            dic_coverage["length"].append(len(gene.sequence))
            pdbs = ""
            for pdbcod, pdb in gene.pdbs.items():
                pdbs += pdbcod + " "

            dic_coverage["pdbs"].append(pdbs)
        pdbs_df = pd.DataFrame.from_dict(dic_coverage)
        return pdbs_df

    def getPdbTaskList(self):
        dic_coverage = {}        
        dic_coverage["gene"] = []
        dic_coverage["accession"] = []                
        dic_coverage["pdb"] = []

        for pdbcod,pdb in self.pdbs.items():            
            dic_coverage["gene"].append(self.gene)
            dic_coverage["accession"].append(self.accession)
            dic_coverage["pdb"].append(pdbcod)                        
        pdbs_df = pd.DataFrame.from_dict(dic_coverage)
        return pdbs_df

    def createPdbConfigYaml(self, pdb, file_path):
        """
        id='1tst'
        chain:'A'
        variant:'Alpha'
        variantfile:'1tst_vars'
        split:'2'
        #???id=3@array=2
        #???id=6@array=3
        repairs:'2'
        #Batch inputs
        jobs:'1234567'
        name:'1tst_2'
        mutation:'.'
        """
        dic_config = {}
        dic_config["id"] = pdb.pdbcode
        dic_config["chain"] = pdb.chain
        dic_config["variant"] = "ANY"
        dic_config["split"] = 2
        dic_config["repairs"] = 2
        dic_config["jobs"] = "1234567"
        dic_config["mutation"] = "."

        import yaml

        with open(file_path, "w") as file:
            outputs = yaml.dump(dic_config, file)
