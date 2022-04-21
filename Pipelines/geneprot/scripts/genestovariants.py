"""
RSA 11.4.22
A few different helper functions to get protein structures from gene names
Via web scraping at different levels of robustness

The main web service we use is bioservices
https://pypi.org/project/bioservices/
https://stackoverflow.com/questions/23443505/how-to-query-uniprot-org-to-get-all-uniprot-ids-for-a-given-species#:~:text=If%20you%20do%20not%20want%20to%20use%20a,%3D%20u.search%20%28%22organism%3A10090%2Band%2Breviewed%3Ayes%22%2C%20columns%3D%22id%2Centry%20name%22%2C%20limit%3D2%29%20print%20%28results%29
"""

import os
import pandas as pd
     
def extractVariantsFromFile(file):
    d = pd.read_csv(file, sep="\t")   
    dic_new_df = {}     
    for i in range(len(d.index)):
    #for i in range(20):
        gene = d["gene"][i].upper()
        
        if gene not in dic_new_df:
            dic_new_df[gene] = {}            
            dic_new_df[gene]['bases'] = []
            dic_new_df[gene]['prev_aa'] = []
            dic_new_df[gene]['residue'] = []
            dic_new_df[gene]['new_aa'] = []
            dic_new_df[gene]['variant'] = []   
               
        variant = d["protein_desc"][i]
        variant_type = d["summary"][i]
        bases = d["mut"][i]    
        if variant_type == "missense":
            if "?" not in variant and "*" not in variant and "delins" not in variant: #? not in exons, * is a stop codon, delins is 2 residues
                if variant not in dic_new_df[gene]['variant']:                    
                    dic_new_df[gene]['bases'].append(len(bases))                        
                    variant = variant[2:]
                    dic_new_df[gene]['variant'].append(variant)
                    aa1 = variant[:1]
                    aa2 = variant[-1:]
                    rid = variant[1:-1]
                    dic_new_df[gene]['prev_aa'].append(aa1)
                    dic_new_df[gene]['residue'].append(int(rid))
                    dic_new_df[gene]['new_aa'].append(aa2)                    
    return dic_new_df
    


    