from Bio import Entrez
import tkinter as tk
import re
import certifi
import urllib3
from bs4 import BeautifulSoup

"""
Created by Andrew Chabot
GUI based webscraping for gene name and minor gene ID analysis built in  
"""
def handleURL(geneURL, geneNum, fileName):
    """
    Takes in a url, turns it to html, and uses helper functions to get a result
    for the most likely X(geneNum) genes in an article
    :param geneURL: URL for the article to search in
    :param geneNum: How many genes to search for in article
    :param fileName: name of article to output info to
    :return: List of gene names
    """
    fle = open(fileName, "w+")
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
    htmlPart = http.request('GET', geneURL)
    print("Finding gene in webpage...")
    fle.write("Finding gene in webpage...\n")
    preArray = createKeyWordPreArr()
    postArray = createKeyWordPostArr()
    soup = BeautifulSoup(htmlPart.data, "html.parser")
    potNames = getPotNames(preArray, 0, soup, geneNum)
    potNames.extend(getPotNames(postArray, 1, soup, geneNum))
    reducedNames = elimNonNames(potNames)
    if len(reducedNames) > 0:
        geneNameList = getMostFreqList(reducedNames, geneNum)
        print("Top Gene: " + geneNameList[0])
        fle.write("Top Gene: " + geneNameList[0] + "\n")
        if(len(geneNameList) < geneNum):
            print("You wanted to find " + geneNum + " genes in this article and I only found " + str(len(geneNameList)))
            fle.write("You wanted to find " + geneNum + " genes in this article and I only found " + str(len(geneNameList)) + "\n")
        if(len(geneNameList) > 1):
            print("Other Genes: ", end="")
            fle.write("Other Genes: ")
        for i in range(len(geneNameList) - 1):
            print(geneNameList[i + 1])
            if i < len(geneNameList) - 2:
                fle.write(", ")
                print(", ", end="")
        print("")
        fle.close()
        return geneNameList
    else:
        if len(potNames) > 0:
            print("Error: Article does not appear to have gene name in text. Might be a video or a general article")
            fle.write("Error: Article does not appear to have gene name in text. Might be a video or a general article\n")
        else:
            print("Error: Does not appear to be an article on genetics")
            print("Error: Does not appear to be an article on genetics\n")
    fle.close()
    return ""

def createKeyWordPreArr():
    """
    creates and returns an array of words that would follow a gene name, will
    be used to search in article's html
    :return: an array of words that would likely follow a gene name
    """
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
    """
    creates and returns an array of words that would be before a gene name,
    will be used to search in article's html
    :return: an array of words that would likely have a gene name follow them
    """
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
    """
    function that takes in an array of possible names and removes them if they
    are in the common dictionary or if they contain any special characters that
    would not tend to be in a gene name
    :param arr: input array of words that will have non potential gene names
    eliminated
    :return: an array of potential gene names
    """
    newArr = []
    dictFile = open('dict.txt', 'r')
    textWall = dictFile.read()
    dictFile.close()
    for ele in arr:
        if not ele.lower().strip(',') in textWall and re.match("^[a-zA-Z0-9,-]*$", ele) and not ele.isdigit():
            newArr.append(ele.strip(','))
    return newArr

def getPotNames(arr, modEq, soup, geneNum):
    """
    function that finds all the potential gene names in an article based on an
    input array
    :param arr: input array of either post or pre arrays that will
    :param modEq: 0 or 1 depending on whether the arr is a post or pre array
    :param soup: html from article URL
    :param geneNum: number of gene names to search for in html
    :return: array of potential names, all that follow or come ahead of the
    given array
    """
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
                            if(len(splat) > i):
                                retArr.append(splat[i])
                else:
                    if i % 2 == modEq:
                        splat = splitHTML[i - 1].split()
                        for i in range(geneNum):
                            if(len(splat) > i):
                                retArr.append(splat[len(splat) - (i + 1)])
    return retArr

def getMostFreqList(arr, geneNum):
    """
    Finds the most frequent names present in the reduced array
    :param arr: reduced array of potential gene names
    :param geneNum: number of gene names being returned
    :return: list of geneNum size that is the most frequent names in the given
    arr
    """
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
    """
    basic function that runs the handleURL function and takes the output and
    puts it into analyzeGene. Uses the text given in the textboxes by the user
    :return: none, prints
    """
    if(e3.get() == ""):
        fileName = "output.txt"
    else:
        fileName = e3.get()
        if fileName[len(fileName) - 4:] != ".txt":
            fileName += ".txt"
    geneNameList = handleURL(e1.get(), int(e2.get()), fileName)
    if(len(geneNameList) > 0):
        analyzeGene(geneNameList, fileName)
    else:
        fle = open(fileName, "a+")
        fle.write("No analysis done, no genes found!")
        print("No analysis done, no genes found!")
        fle.close()

def analyzeGene(geneList, fileName):
    """
    Conducts analysis of all given genes, pritns and outputs information to a
    file
    :param geneList: list of genes to be analyzed
    :param fileName: file to output information to
    :return: none
    """
    fle = open(fileName, "a+")
    email = "afc1755@rit.edu"
    Entrez.email = email
    for i in range(len(geneList)):
        fle.write("Searching for gene: " + geneList[i] + " in NCBI gene database\n")
        print("Searching for gene: " + geneList[i] + " in NCBI gene database")
        handle = Entrez.esearch(db="gene", retmax=5, term=geneList[i])
        dicEle = Entrez.read(handle)
        ids = dicEle.get('IdList')
        printMe = "Top 5 Matching Gene IDs: "
        printMe += str(ids)
        fle.write(printMe + "\n")
        print(printMe)
        for ele in ids:
            fle.write("Info for Gene ID: " + ele + "\n")
            print("Info for Gene ID: " + ele)
            summary = Entrez.esummary(db="gene",id=ele)
            sumRead = Entrez.read(summary)
            handle.close()
            fle.write("Organism: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['Organism']['ScientificName'] + ", " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['Organism']['CommonName'] + "\n")
            fle.write("Name: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['Name'] + "\n")
            fle.write("Desription: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['Description'])
            print("Organism: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['Organism']['ScientificName'] + ", " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['Organism']['CommonName'])
            print("Name: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['Name'])
            print("Desription: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['Description'])
            if(len(sumRead['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo']) > 0):
                fle.write("Chromosome Location: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrLoc'] + "\n")
                fle.write("Chromosome Start and Stop: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStart'] + ", " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStop'] + "\n")
                fle.write("Exon Count: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ExonCount'] + "\n")
                print("Chromosome Location: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrLoc'])
                print("Chromosome Start and Stop: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStart'] + ", " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStop'])
                print("Exon Count: " + sumRead['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ExonCount'])
            else:
                print("Most genetic info unknown about this gene")
                fle.write("Most genetic info unknown about this gene\n")
            fle.write("\n")
            print("\n")
        print("")
        fle.write("\n")
        handle.close()
    fle.close()
    print("Analysis complete!")

#This is all setup of the GUI using TKinter
#I dont know how to do this in its own function, so it is run by default
#all functions are run by clicking on search, a Button that calls analyze
root = tk.Tk()
root.title("Webpage Gene Finder")
label1 = tk.Label(root, fg="maroon", text='URL:')
e1 = tk.Entry(root, width=80)
label2 = tk.Label(root, fg="maroon", text='Number of Genes to Search For:')
e2 = tk.Entry(root, width=20)
label3 = tk.Label(root, fg="maroon", text='File To Save to(Default: here):')
e3 = tk.Entry(root, width=40)
label1.pack()
e1.pack()
label2.pack()
e2.pack()
label3.pack()
e3.pack()
search = tk.Button(root, text='Analyze Article Gene', width=25, command=analyze)
search.pack(side=tk.BOTTOM)
root.mainloop()
