"""
RSA 11.4.22
A few different helper functions to get protein structures from gene names
Via web scraping at different levels of robustness

The main web service we use is bioservices
https://pypi.org/project/bioservices/
https://stackoverflow.com/questions/23443505/how-to-query-uniprot-org-to-get-all-uniprot-ids-for-a-given-species#:~:text=If%20you%20do%20not%20want%20to%20use%20a,%3D%20u.search%20%28%22organism%3A10090%2Band%2Breviewed%3Ayes%22%2C%20columns%3D%22id%2Centry%20name%22%2C%20limit%3D2%29%20print%20%28results%29
"""

import os
from urllib.error import URLError
import requests

from urllib.request import urlretrieve
import Bio.PDB as bio
## Recreating all functions wiuth uniprot API
#https://www.uniprot.org/help/api_queries
def accession_from_uniprot(genename,organism_id,reviewed):
    accessions = []
    url = f"https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:{organism_id}+AND+gene_exact:{genename}"        
    if not reviewed:
        url = f"https://rest.uniprot.org/uniprotkb/search?query=organism_id:{organism_id}+and+reviewed:false+and+gene_exact:{genename}"    
    print(url)      
    # sending get request and saving the response as response object
    ra = requests.get(url=url)    
    # extracting data in json format
    data = ra.json()
    if len(data["results"])> 0:
        for dt in data["results"]:                
            accession = dt["primaryAccession"]      
            accessions.append(accession)
    
    return accessions

def sequence_from_uniprot(accession):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=accession:{accession}"        
    print(url)      
    r = requests.get(url=url)    
    data = r.json()  
    res = data["results"][0]    
    seq = res["sequence"]["value"]    
    return seq


def pdbs_from_accession_uniprot(accession):
    pdb_paths = []  # a tuple of pdb code and path for download
    # first alphafold
    af_model, af_version = "F1", "v2"
    af_name, af_path = getAlphaFoldLink(accession, af_model, af_version)
    # next uniprot PBS experiemntal structures
    url = f"https://rest.uniprot.org/uniprotkb/search?query=accession:{accession}"            
    print(url)      
    pdb_paths.append({"pdb": af_name, "path": af_path})                        
    r = requests.get(url=url)    
    data = r.json()  
    res = data["results"][0]["uniProtKBCrossReferences"]
    for x in res:
        db = x["database"]
        pdb = x["id"]
        if db == "PDB":
            pdb_path = getPDBLink(pdb)
            pdb_paths.append({"pdb": pdb.lower(), "path": pdb_path})                        
    return pdb_paths
    
# ----------------- -------------------------------------------
###       PDB Links                                       ###
# ------------------------------------------------------------
def getPDBeLink(pdb):
    return "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb" + pdb.lower() + ".ent"


def getPDBLink(pdb):
    return "https://files.rcsb.org/download/" + pdb.upper() + ".pdb"


def getAlphaFoldLinkFromPdb(pdbcode):
    ver = pdbcode[-2:]
    mainbit = pdbcode[4:-3]
    bits = mainbit.split("-")
    acc = bits[1]
    model = bits[2]
    return getAlphaFoldLink(acc,model,ver)[1]

def getAlphaFoldLink(accession, model, version):
    af_path = (
        "https://alphafold.ebi.ac.uk/files/AF-"
        + accession.upper()
        + "-"
        + model.upper()
        + "-model_"
        + version.lower()
        + ".pdb"
    )
    af_name = (
        "AF-" + accession.upper() + "-" + model.upper() + "-model_" + version.lower()
    )
    return af_name, af_path


def retrievePdbStructure(url, pdb, path_name):
    if retrieveFile(url, path_name):
        try:
            parser = bio.PDBParser()
            #print("PDB:", path_name)
            struc = parser.get_structure(pdb, path_name)
            return struc
        except:
            if retrieveFile(url, path_name, overwrite=True):
                parser = bio.PDBParser()
                print("PDB:", path_name)
                struc = parser.get_structure(pdb, path_name)
                return struc
    return None


def removePdbStructure(url, pdb, path_name):
    pdb_path = path_name + pdb + ".pdb"
    if os.path.exists(pdb_path):
        os.remove(pdb_path)
    # if os.path.exists(path_name):
    #    os.remove(path_name)


def retrieveFile(url, path_name, overwrite=False):
    if not os.path.exists(path_name) or overwrite:
        try:
            urlretrieve(url, path_name)
            return True
        except URLError as e:
            print("...!!! No data for", url,e)
            return False
    return True


"""
##################### BIOSERVICES IS TOO UNRELIABLE #################################
# ----------------- -------------------------------------------
###       bioservices, uniprot                            ###
# ------------------------------------------------------------
from bioservices import UniProt
def xaccession_from_bioservices(genename,organism_id,reviewed):    
    u = UniProt()

    ## If you have any uniprot problems this line should work so check it, eg you might have a connction error
    #print("QUERY THAT WORKS")
    #res = u.search("P43403", frmt="txt")
    #print(res)
    # You can browse the reults in a browser:
    # https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606

    #search_string = "organism_id:10090+and+reviewed:true+and+gene:" + genename    
    review_string = "yes"
    if not reviewed:
        review_string = "no"

    search_string = "organism:"+organism_id+"+and+reviewed:"+review_string+"+and+gene:" + genename
    print("Searching:", "https://rest.uniprot.org/uniprotkb/search?query=" + search_string)
    result = u.search(
        #"organism:9606+and+reviewed:yes+and+gene:" + genename,
        search_string,
        columns="id,genes",
        limit=5,
    )
    print(result)
    try:
        rows = result.split("\n")
    except:
        raise TypeError("!!! Error in BioServices HTML request !!! Aborting :-(")
    accs = []
    if len(rows) > 1:
        for row in rows:
            print("Row:", row)
            acc_gene = row.split("\t")
            if len(acc_gene) > 1:
                acc = acc_gene[0]
                gn = acc_gene[1].upper()
                gns = gn.split(" ")
                for gn in gns:
                    gn = gn.strip()                
                    if gn == genename:
                        accs.append(acc)
        return accs
    else:        
        return []


def xsequence_from_bioservices(accession):
    u = UniProt()
    seq = u.retrieve(accession, "fasta")
    return seq


def xpdbs_from_accession_bioservices(accession):
    u = UniProt()
    af_model, af_version = "F1", "v2"
    pdb_paths = []  # a tuple of pdb code and path for download
    res = u.mapping("ACC", "PDB_ID", accession)
    af_name, af_path = getAlphaFoldLink(accession, af_model, af_version)
    pdb_paths.append({"pdb": af_name, "path": af_path})
    if accession in res:
        for pdb in res[accession]:
            pdb_path = getPDBLink(pdb)
            pdb_paths.append({"pdb": pdb.lower(), "path": pdb_path})
    return pdb_paths


"""