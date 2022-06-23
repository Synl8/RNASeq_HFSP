from collections import defaultdict
from numpy import mean, log10, std, linspace
import re
from multiprocessing import Pool

from Bio import pairwise2

ggg_test = lambda seq: seq[::-1].startswith("GGGG"[::-1])

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
    seq_x = seq_x_.replace("-", "")
    if test_if_sub(seq_x):
        try:
            alig = pairwise2.align.globalms(ref_seq, seq_x, 2, -2, -4, -1, penalize_end_gaps=False)[0]
        except: 
            return False
        al_pos = "".join([x_n for ref_n, x_n in zip(alig[0], alig[1]) if ref_n != "-"])
        return al_pos[-20:].count("-")/20. < 0.5
    else:
        return False


all_designs = {n.split()[0]: seq for n, seq in read_fasta("./db/pool_1/pool1.fa").items()}
sel_fa_nsu = {n.split()[0]: s for n, s in read_fasta("./blast_res/L422T13_2_p1_cl.fa").items()}

sel_nsu, sel_sub = defaultdict(lambda: []), defaultdict(lambda: [])
# save reads from no-sub condition
for name, seq in read_fasta("./blast_res/L422T13_2_p1_cl.fa").items():
    name_des = name.split("|")[1]
    sel_nsu[name_des] += [seq]

# save active reads
sub_seq_list = list(read_fasta("./blast_res/L422T14_2_p1_cl.fa").items())
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
for (name, seq), ref_seq, act in zip(sub_seq_list, ref_seq_list, res_act):
    name_des = name.split("|")[1]
    if act:
        if name_des in activity_des:
            activity_des[name_des] += 1
        else:
            activity_des[name_des] = 1

for name_des, nb_hits in activity_des.items():
    print(name_des, nb_hits)
