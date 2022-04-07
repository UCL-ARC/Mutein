"""
Getting protein codes from genes
Why does it say a different name
"""
from lxml import html
import requests

page = requests.get('http://econpy.pythonanywhere.com/ex/001.html')
tree = html.fromstring(page.content)
#This will create a list of buyers:
buyers = tree.xpath('//div[@title="buyer-name"]/text()')
#This will create a list of prices
prices = tree.xpath('//span[@class="item-price"]/text()')
print('Buyers: ', buyers)
print('Prices: ', prices)


page = requests.get('https://alphafold.ebi.ac.uk/search/text/NOTCH1?organismScientificName=Homo%20sapiens')

from Bio import Entrez
Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
handle = open("Homo_sapiens.xml", "b")
records = Entrez.parse(handle)
for record in records:
    status = record["Entrezgene_track-info"]["Gene-track"]["Gene-track_status"]
    if status.attributes["value"]=="discontinued":
        continue
    geneid = record["Entrezgene_track-info"]["Gene-track"]["Gene-track_geneid"]
    genename = record["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]
    print(geneid, genename)

