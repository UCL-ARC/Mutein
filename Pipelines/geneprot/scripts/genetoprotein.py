"""
RSA 11.4.22
A few different helper functions to get protein structures from gene names
Via web scraping at different levels of robustness

The main web service we use is bioservices
https://pypi.org/project/bioservices/
https://stackoverflow.com/questions/23443505/how-to-query-uniprot-org-to-get-all-uniprot-ids-for-a-given-species#:~:text=If%20you%20do%20not%20want%20to%20use%20a,%3D%20u.search%20%28%22organism%3A10090%2Band%2Breviewed%3Ayes%22%2C%20columns%3D%22id%2Centry%20name%22%2C%20limit%3D2%29%20print%20%28results%29
"""

from bs4 import BeautifulSoup  # https://beautiful-soup-4.readthedocs.io/en/latest/#quick-start
import requests
from bioservices import UniProt

      
#------------------------------------------------------------
###       bioservices, uniprot                            ###
#------------------------------------------------------------
def accession_from_bioservices(genename):
    
    u = UniProt()    
    result = u.search("organism:9606+and+reviewed:yes+and+gene:" + genename, columns="id", limit=1)
    rows = result.split("\n")    
    if len(rows)>1:
        return rows[1] #don't return the header
    else:
        print("(!)" +result)
        return ""
def pdbs_from_accession_bioservices(accession):
    u = UniProt()
    res = u.mapping("ACC", "PDB_ID", accession)        
    return res[accession]
#------------------------------------------------------------
###       bioservices, uniprot                            ###
#------------------------------------------------------------
        

    