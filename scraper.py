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
    preArray = createKeyWordPreArr()
    postArray = createKeyWordPostArr()
    soup = BeautifulSoup(htmlPart.data, "html.parser")
    potNames = getPotNames(preArray, 0, soup)
    potNames.extend(getPotNames(postArray, 1, soup))
    reducedNames = elimNonNames(potNames)
    if len(reducedNames) > 0:
        geneName = getMostFreq(reducedNames)
        print("Best Guess for Gene: " + geneName)
    else:
        if len(potNames) > 0:
            print("Article does not appear to have gene name in text. Might be a video or a general article")
        else:
            print("Does not appear to be an article on genetics")

def createKeyWordPreArr():
    kwArr = []
    kwArr.append(" gene ")
    kwArr.append(" mutations ")
    kwArr.append(" mutation ")
    kwArr.append(" expression ")
    return kwArr

def createKeyWordPostArr():
    kwArr = []
    kwArr.append(" gene called ")
    kwArr.append(" mutation in ")
    kwArr.append(" mutation in the ")
    kwArr.append(" variants of ")
    kwArr.append(" variance in ")
    kwArr.append(" named ")
    kwArr.append(" expression of ")
    kwArr.append(" overexpression of ")
    kwArr.append(" underexpression of ")
    return kwArr

def elimNonNames(arr):
    newArr = []
    dictFile = open('dict.txt', 'r')
    textWall = dictFile.read()
    dictFile.close();
    for ele in arr:
        if not ele in textWall:
            newArr.append(ele)
    return newArr

def getPotNames(arr, modEq, soup):
    retArr = []
    for term in arr:
        if term in str(soup):
            splitHTML = str(soup).split(term)
            for i in range(0, len(splitHTML)):
                if modEq == 1:
                    if i % 2 == modEq:
                        retArr.append(splitHTML[i].split()[0])
                else:
                    if i % 2 == modEq:
                        splitAgain = splitHTML[i].split()
                        retArr.append(splitAgain[len(splitAgain) - 1])
    return retArr

def getMostFreq(arr):
    freqDict = {}
    for term in arr:
        if term in freqDict:
            freqDict[term] = freqDict[term] + 1
        else:
            freqDict[term] = 1
    return max(freqDict)

def main():
    geneURL = input("Please enter a url: ")
    handleURL(geneURL)

main()