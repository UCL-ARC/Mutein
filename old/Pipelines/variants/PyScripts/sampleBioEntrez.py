#https://bioinformatics.stackexchange.com/questions/19356/downloading-genomic-protein-files-from-accessions-in-python
"""
from Bio import Entrez
Entrez.email = "m@M__"
#handle = Entrez.efetch(db="protein", id='"GCF_000145295", "GCF_000145294", "GCF_000145293"', retmode="text",rettype="gb")
handle = Entrez.efetch(db="protein")
records = Entrez.parse(handle)
handle.close()
for record in records:
    print (record)
"""

#https://biopython.org/docs/1.76/api/Bio.Entrez.html
#https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
from Bio import Entrez
Entrez.email = "Your.Name.Here@example.org"
#handle = Entrez.einfo() # or esearch, efetch, ...
#handle = Entrez.esummary(db="pubmed", id="19304878,14630660", retmode="xml")
#handle = Entrez.efetch(db="protein", id="15718680,NP_001098858.1,119703751", retmode="xml")
handle = Entrez.efetch(db="nuccore", id="21614549", strand="1", seq_start="1", seq_stop="100", retmode="xml")
#handle = Entrez.efetch(db="protein", id="P46531", retmode="xml")
record = Entrez.read(handle)
handle.close()
print("#################################")
count = 0
for r in record:
    count +=1
    for rr,pp in r.items():
        print(count,rr)
    print(r["GBSeq_locus"])
    print(r["GBSeq_length"])
    print(r["GBSeq_sequence"])
    