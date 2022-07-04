"""
RSA 30/6/22
------------------------
Class to manage gene and chromosome position

"""
import os
import requests

class Chme:
    def __init__(self, chme,oranism_id):
        self.chme = chme
        self.org_id = oranism_id
        self.genes = {}
        self.gene_seqs = {}
                        
    def addGene(self, gene,start,end):
        if gene in self.genes:            
            raise Exception("Can't have genes in multiple regions")        
        self.genes[gene] = [start,end]

        
    def findGene(self,pos):
        for gene,startend in self.genes.items():            
            if pos >= startend[0] and pos <= startend[1]:
                return gene
        return ""

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
        
