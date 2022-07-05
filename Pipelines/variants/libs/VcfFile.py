"""
RSA 5/7/22
------------------------
Class to manage vcf file

"""
import allel
import pandas as pd

class VcfFile:
    def __init__(self, vcf_file):
        self.vcf_file = vcf_file                

    def getVariants(self):
        dic_vcf = {}
        callset = allel.read_vcf(self.vcf_file,fields="*")        
        vshp = callset["variants/POS"].shape
        print(vshp)        
        for key,arr in callset.items():                                    
            shp = arr.shape
            if "variants" in key:
                key = key[9:]
            if shp == vshp:                
                dic_vcf[key] = arr                            
        # make it a dataframe
        df_vcf = pd.DataFrame.from_dict(dic_vcf)   
        #add the mutations too - as it is snp we can take the first column
        alts = callset["variants/ALT"]
        muts = []
        for alt in alts:
            muts.append(alt[0])
        df_vcf["ALT"] = muts
        df_vcf = df_vcf.query("`is_snp` == True")
        df_vcf_cols = df_vcf[["CHROM","REF","ALT","POS"]]
        return df_vcf_cols  
        