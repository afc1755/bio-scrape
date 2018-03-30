import Bio
from Bio import Entrez
from Bio import SeqIO
import certifi
import urllib3
from bs4 import BeautifulSoup

"""Created by Andrew Chabot
command line webscraping, main center for project code
"""

def handleURL(geneURL):
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    htmlPart = http.request('GET', geneURL)
    keywordArray = [];
    potentialNamesArr = [];
    soup = BeautifulSoup(htmlPart.data, "html.parser")
    """if "gene called" in str(soup):
        splitHTML = str(soup).split("gene called")
        potGeneName = splitHTML[1].split()[0]
    """
    #if "gene is called" in soup:

    #print("Best Guess for Gene beign talked about:" + geneName)

def main():
    geneURL = input("Please enter a url: ")
    handleURL(geneURL)

main()