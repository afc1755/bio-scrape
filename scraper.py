import Bio
from Bio import Entrez
from Bio import SeqIO
import certifi
import urllib3
from bs4 import BeautifulSoup

"""Created by Andrew Chabot
command line webscraping, main center for project code
"""

def main():
    #geneName = input("Please enter a url")
    page = 'http://www.bloomberg.com/quote/INDU:IND'
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    htmlPart = http.request('GET', page)
    #print(htmlPart)
    soup = BeautifulSoup(htmlPart.data, "html.parser")
    #print(soup)
    name_box = soup.find('h1', attrs = {'class': 'name'})
    #print(name_box)
    name = name_box.text.strip()
    print(name)
    price_box = soup.find('div', attrs = {'class': 'price'})
    price = price_box.text
    print(price)
    """Entrez.email = "afc1755@rit.edu"
    with Entrez.efetch(db="gene", rettype="fasta", retmode="text", id="6273291") as handle:
        seq_record = SeqIO.read(handle, "fasta")
    print("%s with %i features" % (seq_record.id, len(seq_record.features)))
    """
main()