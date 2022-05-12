"""
RSA 12/5/22
------------------------
Class to manage coverage of a pdb and the data set it needs

"""
import os
from shutil import copyfile
import pandas as pd


class PdbRunner:
    def __init__(self,pdbcode):
        self.pdbcode = pdbcode
    
    def copyToInput(self, thruput_dir, input_dir, df):
        # first copy the pdb file to the pdb directoru
        thru_file = thruput_dir + self.pdbcode.lower() + ".pdb"
        in_file = input_dir + self.pdbcode.lower() + ".pdb"
        cov_file = input_dir + "Coverage.csv"
        print("### copying file", thru_file,in_file)
        copyfile(thru_file,in_file)

        # Then create a pdb dataframe for just that one pdb
        df_one = df.query("pdb == '" + self.pdbcode + "'")
        df_one.to_csv(cov_file,index=False)

    def getVariantCandidatesDataFrame(self,gene,accession,variants, cov_df):
        dic_variants = {}
        dic_variants["gene"] = []
        dic_variants["accession"] = []
        dic_variants["variant"] = []
        dic_variants["residue"] = []
        dic_variants["bases"] = []
        dic_variants["candidates"] = []

        for vrcod, vr in variants.items():
            candidate = ""
            # print(vr.variant)            
            if self.matchesResidue(vr.residue,cov_df):
                candidate += self.pdbcode + " "
            dic_variants["residue"].append(vr.residue)
            dic_variants["gene"].append(gene)
            dic_variants["accession"].append(accession)
            dic_variants["candidates"].append(candidate)
            dic_variants["variant"].append(vr.variant)
            dic_variants["bases"].append(vr.bases)

        pdbs_df = pd.DataFrame.from_dict(dic_variants)
        pdbs_df = pdbs_df.sort_values(by="residue", ascending=True)
        return pdbs_df

        