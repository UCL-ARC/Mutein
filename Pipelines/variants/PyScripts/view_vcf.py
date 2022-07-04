"""
RSA 30/6/22
Use the packakge scikit-allel
(pip install scikit-allel)
https://scikit-allel.readthedocs.io/en/stable/
Some instructions here:
http://alimanfoo.github.io/2017/06/14/read-vcf.html
"""

import sys
import os


dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]
lib_path = "/".join(dirs) + "/libs"
sys.path.append(lib_path)
import Chme


import allel
import pandas as pd

## INPUTS ##
hardcoded_vcf_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/keogh2018.vcf"
hardcoded_gene_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/gene_list"
hardcoded_interval_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/intervals.list"
hardcoded_organism_id = "/home/rachel/UCL/github/MuteinData/vcf_keogh/organism_id.txt"
## outputs ## 
hardcoded_out_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/info.csv"
print("################ PARSE VCF FILE ###############")
with open(hardcoded_vcf_path, mode='r') as org:
    print(org.readlines()[:50])

dic_vcf = {}
callset = allel.read_vcf(hardcoded_vcf_path,fields="*")
length = 0
for key,arr in callset.items():                
    shp = arr.shape
    if len(shp) > 1:
        print(shp)        
    else:
        if length == 0:
            length = shp
            print(shp)
        print(arr)
        if shp == length:
            dic_vcf[key] = arr
        else:
            print(key, shp)

# make it a dataframe
df_vcf = pd.DataFrame.from_dict(dic_vcf)   
print("Save VCF dataframe to",hardcoded_out_path)
df_vcf.to_csv(hardcoded_out_path,index=False)


