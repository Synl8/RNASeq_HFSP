settings = {
"databaseFilePath" : "pool0_1.fa",
"subFilePath" : "L422T14_cl.fa",
"cleanSubFilePath" : "L422T14_cl.fa",
# "cleanSubFilePath" : "L422T14_cl_cleaned.fa",
"nSubFilePath" : "L422T13_cl.fa",
"cleanNSubFilePath" : "L422T13_cl.fa",
# "cleanNSubFilePath" : "L422T13_cl_cleaned.fa",
"mafftPath" : "/usr/bin/mafft",
"exon" : "aatccgttggtgctg",
"sub2": "AACTTCAAATATCTTCGGAACTCA",
"3primeCleanerIdentifier" : "GATCGGAAGAGCACA",
"5primeCleanerIdentifier" : "gtgactggagttcagacgtgtgctcttccgatc",
"G_test" : "[G]{5,10}$",
}


def getSettings():
    return settings 