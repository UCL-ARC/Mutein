"""
RSA 30/6/22
------------------------
Class to manage gene and chromosome position

"""
import os
import requests

class UniProt:
    def __init__(self, oranism_id):
        #self.chme = chme
        self.org_id = oranism_id
        self.genes = {}
        self.gene_seqs = {}
                            
    def seqFromUniProt(self,gene):    
        if gene not in self.gene_seqs:
            url = f"https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:{self.org_id}+AND+gene_exact:{gene}"                        
            print(url)              
            ra = requests.get(url=url)            
            data = ra.json()        
            if len(data)== 0:
                url = f"https://rest.uniprot.org/uniprotkb/search?query=reviewed:false+AND+organism_id:{self.org_id}+AND+gene_exact:{gene}"        
                print(url)              
                ra = requests.get(url=url)            
                data = ra.json()        
            elif len(data["results"])== 0:
                url = f"https://rest.uniprot.org/uniprotkb/search?query=reviewed:false+AND+organism_id:{self.org_id}+AND+gene_exact:{gene}"        
                print(url)              
                ra = requests.get(url=url)            
                data = ra.json()       
            seq,length,accessions = "",0,[]
            if len(data["results"])> 0:
                for dt in data["results"]:                
                    accession = dt["primaryAccession"]      
                    accessions.append(accession)                        
                res = data["results"][0]    
                seq = res["sequence"]["value"]
                length = res["sequence"]["length"]
                    
            self.gene_seqs[gene] = seq
        else:
            seq = self.gene_seqs[gene]
        return seq
        
