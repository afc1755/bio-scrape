from Bio import Entrez
import tkinter as tk
import re
import certifi
import urllib3
from bs4 import BeautifulSoup

"""Created by Andrew Chabot
command line webscraping, main center for project code
"""

def handleURL(geneURL, geneNum):
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    htmlPart = http.request('GET', geneURL)
    preArray = createKeyWordPreArr()
    postArray = createKeyWordPostArr()
    soup = BeautifulSoup(htmlPart.data, "html.parser")
    potNames = getPotNames(preArray, 0, soup, geneNum)
    potNames.extend(getPotNames(postArray, 1, soup, geneNum))
    reducedNames = elimNonNames(potNames)
    if len(reducedNames) > 0:
        geneNameList = getMostFreqList(reducedNames, geneNum)
        print("Top Gene: " + geneNameList[0])
        if(len(geneNameList) < geneNum):
            print("You wanted to find " + geneNum + " genes in this article and I only found " + len(geneNameList))
        if(len(geneNameList) > 1):
            print("Other Genes: ", end = "")
        for i in range(len(geneNameList) - 1):
            print(geneNameList[i + 1])
            if i < len(geneNameList) - 2:
                print(", ", end = "")
        return geneNameList
    else:
        if len(potNames) > 0:
            print("Article does not appear to have gene name in text. Might be a video or a general article")
        else:
            print("Does not appear to be an article on genetics")
    return ""

def createKeyWordPreArr():
    kwArr = []
    kwArr.append(" gene")
    kwArr.append(" is a gene")
    kwArr.append(" mutations")
    kwArr.append(" mutation")
    kwArr.append(" expression")
    kwArr.append(" splicing")
    kwArr.append(", a gene")
    kwArr.append(" effects")
    return kwArr

def createKeyWordPostArr():
    kwArr = []
    kwArr.append("genes called ")
    kwArr.append("gene called ")
    kwArr.append("gene, ")
    kwArr.append("alleles for ")
    kwArr.append("attributed to ")
    kwArr.append("mutation in ")
    kwArr.append("mutations in ")
    kwArr.append("mutation in the ")
    kwArr.append("variants of ")
    kwArr.append("variation of ")
    kwArr.append("variance in ")
    kwArr.append("allele in ")
    kwArr.append("known as ")
    kwArr.append("named ")
    kwArr.append("expression of ")
    kwArr.append("overexpression of ")
    kwArr.append("underexpression of ")
    kwArr.append("gene ")
    kwArr.append("genes ")
    kwArr.append("genes like ")
    kwArr.append(", or ")

    return kwArr

def elimNonNames(arr):
    newArr = []
    dictFile = open('dict.txt', 'r')
    textWall = dictFile.read()
    dictFile.close()
    for ele in arr:
        if not ele.lower().strip(',') in textWall and re.match("^[a-zA-Z0-9,-]*$", ele) and not ele.isdigit():
            newArr.append(ele.strip(','))
    return newArr

def getPotNames(arr, modEq, soup, geneNum):
    retArr = []
    soupStr = str(soup)
    for term in arr:
        if term.upper() in soupStr.upper():
            splitHTML = re.split(term, soupStr, flags=re.IGNORECASE)
            for i in range(0, len(splitHTML)):
                if modEq == 1:
                    if i % 2 == modEq:
                        splat = splitHTML[i].split()
                        for i in range(geneNum):
                            if(len(splat) >= i):
                                retArr.append(splat[i])
                else:
                    if i % 2 == modEq:
                        splat = splitHTML[i - 1].split()
                        for i in range(geneNum):
                            if(len(splat) >= i):
                                retArr.append(splat[len(splat) - (i + 1)])
    return retArr

def getMostFreqList(arr, geneNum):
    freqDict = {}
    for term in arr:
        if term in freqDict:
            freqDict[term] += 1
        else:
            freqDict[term] = 1
        if len(term) == 5:
            freqDict[term] += 3
        if term.upper() == term:
            freqDict[term] += 3
    lst = []
    for i in range(geneNum):
        if(len(freqDict) > 0):
            lst.append(max(freqDict, key=lambda key: freqDict[key]))
            freqDict.pop(max(freqDict, key=lambda key: freqDict[key]))
    return lst

def analyze():
    geneNameList = handleURL(e1.get(), int(e2.get()))
    analyzeGene(geneNameList)

def analyzeGene(geneList):
    email = "afc1755@rit.edu"
    Entrez.email = email
    for i in range(len(geneList)):
        print("Analyzing gene: " + geneList[i])
        handle = Entrez.esearch(db="gene", term=geneList[i])
        printMe = "Start"
        while printMe != "":
            printMe = handle.readline().strip()
            print(printMe)
        print("Analysis complete!")

root = tk.Tk()
root.title("Genetic Article Analysis")
label1 = tk.Label(root, fg="dark green", text='Enter URL:')
e1 = tk.Entry(root, width=80)
label2 = tk.Label(root, fg="dark green", text='Enter Number of Genes to Search For:')
e2 = tk.Entry(root, width=80)
label1.pack()
e1.pack()
label2.pack()
e2.pack()
search = tk.Button(root, text='Analyze Article Gene', width=25, command=analyze)
search.pack(side=tk.BOTTOM)
root.mainloop()
