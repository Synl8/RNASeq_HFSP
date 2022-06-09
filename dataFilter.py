from settings import *
from collections import defaultdict
from numpy import append, mean, log10, std, linspace
import re, sys, os, subprocess, time
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
import json
import pandas as pd
# import tkinter as tk
import numpy  as np


################################################################################
# global variables
ggg_test = lambda seq: seq[::-1].startswith("GGGG"[::-1])
settings = getSettings()

################################################################################

# functions

def saveData(dict, file):
    data_file = open(file, "w")
    json.dump(dict, data_file)
    data_file.close()

def loadData(file):
    with open(file) as json_file:
        data = json.load(json_file)
        return data

def read_fasta(infile):
    "read fasta files, return a dictionary"
    results = {}
    for l in open(infile):
        if l.startswith(">"):
            name = l[1:].strip()
            results[name] = ""
        else:
            results[name] += l.strip()
    return results

def test_if_sub(seq):
    return re.search(".+GGAACTTCAAATATCTTCGG.?", seq) is not None

def test_activity(args):
    ref_seq, seq_x_ = args
    seq_x = seq_x_.replace("-", "").replace("U", "T")
    if test_if_sub(seq_x):
        try:
            seq_x = seq_x[:seq_x.find("GGAACTTCAAATATCTTCGG")]
            alig = pairwise2.align.globalms(ref_seq, seq_x, 2, -2, -4, -1, penalize_end_gaps=False)[0]
        except: 
            return False
        al_pos = "".join([x_n for ref_n, x_n in zip(alig[0], alig[1]) if ref_n != "-"])
        return al_pos[-20:].count("-")/20. < 0.5
    else:
        return False


# usefull global variables 
all_designs = {n.split()[0]: seq for n, seq in read_fasta(settings["databaseFilePath"]).items()}
sel_fa_nsu = {n.split()[0]: s for n, s in read_fasta(settings["cleanSubFilePath"]).items()}
sel_nsu, sel_sub = defaultdict(lambda: []), defaultdict(lambda: [])

def filter():
    # save reads from no-sub condition
    for name, seq in read_fasta(settings["cleanNSubFilePath"]).items():
        name_des = name.split("|")[1]
        sel_nsu[name_des] += [seq]

    # save active reads
    sub_seq_list = list(read_fasta(settings["cleanSubFilePath"]).items())
    ref_seq_list = []

    for name, seq in sub_seq_list:
        name_des = name.split("|")[1]
        ref_des = all_designs[name_des]
        ref_seq_list += [ref_des]

    pool = Pool(100) 

    # test
    res_act = pool.map(test_activity, ((ref_seq, seq) for (name, seq), ref_seq in
                                zip(sub_seq_list, ref_seq_list)))

    activity_des = {}
    nsu_count_des = {}

    for name, seq in sel_nsu.items():
        nsu_count_des[name] = len(seq)


    for (name, seq), ref_seq, act in zip(sub_seq_list, ref_seq_list, res_act):
        name_des = name.split("|")[1]
        if act:
            if name_des in activity_des:
                activity_des[name_des] += 1
            else:
                activity_des[name_des] = 1

            
    ratioList = []
    nbNsuList = []
    nbActifList = []
    nbMutationList = []
    G_testList= []
    nameList = list(activity_des.keys())

    for name in nsu_count_des.keys():
        if name not in nameList : nameList.append(name)

    azo_seq = str.upper(all_designs["AZOARCUS"]).replace("U", "T")
    for name in nameList:
        nbActif =0
        ref_seq = str.upper(all_designs[name]).replace("U", "T")
        nbNsu = 0
        ratio = pd.NA
        nbMut = lambda azo_seq, ref_seq: sum(ei != ej for ei, ej in zip(azo_seq, ref_seq))
        if name in activity_des.keys() : nbActif= activity_des[name]
        if name in nsu_count_des.keys() : nbNsu= nsu_count_des[name]
        if nbActif > 0 and nbNsu > 0 : ratio = nbActif / nbNsu
        nbActifList.append(nbActif)
        nbNsuList.append(nbNsu)
        ratioList.append(ratio)
        nbMutationList.append(nbMut(ref_seq, azo_seq))
        G_testList.append(re.search("[G]{4,8}$", ref_seq) is not None)



    csvData = {"VariantName":nameList, "nsubCountRead":nbNsuList, "subCountRead": nbActifList, "ActRate": ratioList, "Nb_Mut":nbMutationList, "G_test":G_testList}
    df = pd.DataFrame.from_dict(csvData)
    df = df[df["ActRate"] != pd.NA]
    if not os.path.isdir("results_csv") : os.makedirs("results_csv")
    df.to_csv("results_csv/"+ settings["databaseFilePath"].split(".")[0].replace(".", "") + "_activity.csv", index=False,   sep=",")
    saveData({"data":(activity_des, sub_seq_list, ref_seq_list, res_act)}, "save.json")
    return activity_des, sub_seq_list, ref_seq_list, res_act

def printActiveVariant(activity_des):
    for name_des, nb_hits in activity_des.items():
        print(name_des, nb_hits)


def askVariantToCheck(activity_dict, sub_seq_list, ref_seq_list, res_act):
    variantName = input("Please enter the name of the variant you want to check : ")
    while not variantName in activity_dict.keys():
        print("List of available names : ")
        print(activity_dict.keys())
        variantName = input("Please choose a variant in the list above : ")

    return alignVariant(variantName, activity_dict, sub_seq_list, ref_seq_list, res_act)

def alignVariant(variantName, activity_dict, sub_seq_list, ref_seq_list, res_act):
    recordList = []
    for (name, seq), ref_seq, act in zip(sub_seq_list, ref_seq_list, res_act):
        name_des = name.split("|")[1]
        if act and name_des == variantName:
            if len(recordList) == 0: 
                ref_seq += settings["exon"] + settings["sub2"]
                recordList.append(SeqRecord(Seq(ref_seq), id = variantName, description="", name=""))
            recordList.append(SeqRecord(Seq(seq), id = name, description="", name=""))
    
    if not os.path.isdir("results") : os.makedirs("results")
    
    fastaFile ="results/" + variantName+ ".fa"
    alignedFile = "results/" + variantName+ "_aln.fa"
    SeqIO.write(recordList , fastaFile, "fasta")
    mafft_cline = settings["mafftPath"] + " --auto --out "+ os.path.abspath(alignedFile) + " " + os.path.abspath(fastaFile)
    print("Sometimes the subproccess get locked for no reason. If it happens, please quit and run the following command manually :")
    print(mafft_cline)
    print("Processing alignment")
    child= subprocess.Popen(str(mafft_cline), stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32"))
    child.wait()
    
    if os.path.isfile(alignedFile): 
        print("Alignment finished : the aligned file is " + alignedFile)
    else: 
        print("It seems that an error occured, the alignement could not be done.")
    
    if os.path.isfile(fastaFile): 
        os.remove(fastaFile)

    return input("Type exit to exit or enter to check for another variant.")




def main():
    print("Applying filter.py")
    t1 = time.time()
    if not os.path.isfile("save.json"):
        activity_dict, sub_seq_list, ref_seq_list, res_act = filter()
        t2 = time.time()
        # finishPopUpDialogBox("Filtering fnish", "Work done")
    else : 
        data = loadData("save.json")
        activity_dict, sub_seq_list, ref_seq_list, res_act = data["data"]
        t2 = time.time()
    print("Filtering finished in " + str(t2-t1) + " seconds.")
    action = askVariantToCheck(activity_dict, sub_seq_list, ref_seq_list, res_act)
    while action != "exit":
        action = askVariantToCheck(activity_dict, sub_seq_list, ref_seq_list, res_act)


def minEditDistanceV2(s1, s2):
# initialisation  of a matrix full of 0, len()+1, since first character is " ", therefore 0
    matrix = np.zeros ((len(s1)+1, len(s2)+1)) 
	# setting first column to the word, the length is +1 since first character is empty
    matrix [0:len(s1)+1, 0] = [x for x in range (0, len(s1)+1)] 
	# setting first row of matrix to the word, the length is +1 since the first character is empty
    matrix [0,0: len(s2)+1] = [y for y in range (0, len(s2)+1)] 

# iterating over columns in a row, starting from index 1, since  0 was previously preset
    for x in range(1, len(s1)+1): # iterating rows, we start from 1 because 0 was preset before
        for y in range(1, len(s2)+1): 			# if characters are equal, we do nothing
			# checking if the char characters are equal, if yes, nothing happen
            if s1[x-1] == s2[y-1]: 
                matrix [x,y] = matrix[x-1,y-1]
				# if they are different, checking which opperation is the least costly: substitution, deletion or insertion
            else: 
                matrix [x,y] = min(
                    matrix[x-1,y] + 1,
                    matrix[x-1,y-1] + 1,
                    matrix[x,y-1] + 1
                )
	#printing the matrix
    # print (matrix) 

	# returning  the last element of the matrix, corresponding to our minimal edit distance

    return (matrix[len(s1), len(s2)])

if __name__ == '__main__':
    main()