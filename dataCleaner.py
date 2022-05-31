from settings import *
import re, os, time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

settings = getSettings()
def getFastaIntoDict(fasta):
    seqDict = {}
    for record in SeqIO.parse(fasta, "fasta"):
        seqDict[record.id] = str(record.seq)

    return seqDict



def cleanFasta(fastaFile, referenceDict):
    threePrimeCleanerIdentifier = str.upper(settings["3primeCleanerIdentifier"])
    fivePrimeCleanerIdentifier = str.upper(settings["3primeCleanerIdentifier"])
    recordList = []
    nb5PrimeCleaned = 0
    nb3PrimeCleaned = 0
    for record in SeqIO.parse(fastaFile, "fasta"):
        id = record.id
        seq = str(record.seq)
        seq = str.upper(seq)
        if "|" in id : id = id.split("|")[1]

        if threePrimeCleanerIdentifier in seq and not threePrimeCleanerIdentifier in referenceDict[id]:
            seq= seq[seq.find(threePrimeCleanerIdentifier) + len(threePrimeCleanerIdentifier):]
            nb3PrimeCleaned+=1

        if fivePrimeCleanerIdentifier in seq and not fivePrimeCleanerIdentifier in referenceDict[id]: 
            # print(record.seq)
            seq = seq[:seq.find(fivePrimeCleanerIdentifier)]
            nb5PrimeCleaned+=1
            # print("After cut of 5s")
            # print(record.seq)
            # return
        
        record.seq = Seq(seq)
        recordList.append(record)

    SeqIO.write(recordList, fastaFile.replace(".fa", "_cleaned.fa"), "fasta")
    print("in the file :", fastaFile)
    print(str(nb5PrimeCleaned) + " 5prime have been removed")
    print(str(nb3PrimeCleaned) + " 3prime have been removed")

def cleanAlignementFasta(fastaFile, referenceDict):
    threePrimeCleanerIdentifier =str.upper(settings["3primeCleanerIdentifier"])
    fivePrimeCleanerIdentifier =str.upper(settings["3primeCleanerIdentifier"])
    recordList = []
    nb5PrimeCleaned = 0
    nb3PrimeCleaned = 0
    for record in SeqIO.parse(fastaFile, "fasta"):
        record.seq = Seq(str.upper(str(record.seq).replace("-", "")))
        # print(record.seq)
        id = record.id
        if "|" in id : id = id.split("|")[1]
        if fivePrimeCleanerIdentifier in record.seq and not fivePrimeCleanerIdentifier in referenceDict[id]: 
            record.seq = record.seq[str(record.seq).find(fivePrimeCleanerIdentifier)+ len(fivePrimeCleanerIdentifier):]
            nb5PrimeCleaned+=1
        if threePrimeCleanerIdentifier in record.seq and not threePrimeCleanerIdentifier in referenceDict[id]:
            record.seq = record.seq[:str(record.seq).find(threePrimeCleanerIdentifier)]
            nb3PrimeCleaned+=1
        recordList.append(record)

    SeqIO.write(recordList, fastaFile.replace(".fa", "_cleaned.fa"), "fasta")
    print("in the file :", fastaFile)
    print(str(nb5PrimeCleaned) + " 5prime have been removed")
    print(str(nb3PrimeCleaned) + " 3prime have been removed")


def main():
    print("getting reference Dict")
    t1 = time.time()
    referenceDict = getFastaIntoDict(settings["databaseFilePath"])
    t2 = time.time()
    print("refDict was collected in " + str(t2-t1) + " secondes")
    # cleanAlignementFasta("results/DCA_SB_75_460_aln.fa", referenceDict)
    print("start cleaning : ", settings["subFilePath"])
    t1 = time.time()
    cleanFasta(settings["subFilePath"], referenceDict)
    t2 = time.time()
    print("cleaning for " + settings["subFilePath"] + " was made in " + str(t2-t1) + " secondes")
    print("start cleaning : ", settings["nSubFilePath"])
    t1= time.time()
    cleanFasta(settings["nSubFilePath"], referenceDict)
    t2 = time.time()
    print("cleaning for " + settings["nSubFilePath"] + " was made in " + str(t2-t1) + " secondes")
    print("finish")

if __name__ == '__main__':
    main()