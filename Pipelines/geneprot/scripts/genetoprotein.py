"""
RSA 11.4.22
A few different helper functions to get protein structures from gene names
Via web scraping at different levels of robustness

"""


def proteins_from_genes(genename, method):
    if method == "uniprot":
        return proteins_from_uniprot(genename)
    else:
        return []


def proteins_from_uniprot(genename):
    from bs4 import (
        BeautifulSoup,
    )  # https://beautiful-soup-4.readthedocs.io/en/latest/#quick-start

    page = requests.get("https://www.ncbi.nlm.nih.gov/gene/?term=notch1")
    soup = BeautifulSoup(page.text)
    lis = soup.find_all("div", class_="ncbi-doc-id")
    for li in lis:
        print("classa", str(li))

    page = requests.get(
        "https://www.uniprot.org/uniprot/?query=notch1&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score"
    )
    soup = BeautifulSoup(page.text)
    lis = soup.find_all("td", class_="entryID")
    for li in lis:
        print("classb", str(li).split("/"))
