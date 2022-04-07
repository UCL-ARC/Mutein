"""
Getting protein codes from genes
Web scraping
???
"""
from xml.dom.minidom import ElementInfo
from lxml import html
import requests

# This shows how it wokrs on a tutorial site
page = requests.get('http://econpy.pythonanywhere.com/ex/001.html')
tree = html.fromstring(page.content)
#This will create a list of buyers:
buyers = tree.xpath('//div[@title="buyer-name"]/text()') # <div title="buyer-name">Patty Cakes</div>
#This will create a list of prices
prices = tree.xpath('//span[@class="item-price"]/text()') #<span class="item-price">$15.26</span><br>
print('Buyers: ', buyers)
print('Prices: ', prices)

from bs4 import BeautifulSoup
# https://beautiful-soup-4.readthedocs.io/en/latest/#quick-start
page = requests.get('https://www.ncbi.nlm.nih.gov/gene/?term=notch1')
soup = BeautifulSoup(page.text)
lis = soup.find_all('div',class_="ncbi-doc-id")
for li in lis:
    print('classa',str(li))

page = requests.get("https://www.uniprot.org/uniprot/?query=notch1&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score")
soup = BeautifulSoup(page.text)
lis = soup.find_all('td',class_="entryID")
for li in lis:
    print('classb',str(li).split('/'))
  


# unfortunately alphafold is all javascript
page = requests.get('https://alphafold.ebi.ac.uk/search/text/NOTCH1?organismScientificName=Homo%20sapiens')

from Bio import Entrez
Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
"""handle = open("Homo_sapiens.xml", "b")
records = Entrez.parse(handle)
for record in records:
    status = record["Entrezgene_track-info"]["Gene-track"]["Gene-track_status"]
    if status.attributes["value"]=="discontinued":
        continue
    geneid = record["Entrezgene_track-info"]["Gene-track"]["Gene-track_geneid"]
    genename = record["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]
    print(geneid, genename)"""

