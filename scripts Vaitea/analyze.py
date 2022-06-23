from src.seq_io import read_fasta, comp_seq
from glob import glob
from random import uniform
import re
from multiprocessing import Pool

all_res = {}


def read_blast_output(infile):
    results = {}
    for l in open(infile):
        qseqid, sseqid, pident, qcovs, qcovus, strand, sseq, qseq = l.strip().split()
        results[qseqid] = (qseqid, sseqid, float(pident), int(qcovs), int(qcovus), strand, sseq, qseq)
    return results


def check_func(args):
    bl_out = args
    samp_n = re.search("./blast_res\/(.+)\.bl", bl_out).groups()[0]
    bl_dat = read_blast_output(bl_out)
    samp_fa = read_fasta(f"./reads/{samp_n}.fa")
    nb_seq = 0
    out_tab_res = {} 
    with open(f"./blast_res/{samp_n}_cl.fa", "w") as out_fasta, open(f"./blast_res/{samp_n}_cl.tab", "w") as out_tab:
        for qseqid, sseqid, pident, qcovs, qcovus, strand, sseq, qseq in bl_dat.values():
            if qcovs >= 50 and pident == 100:
                if strand == "minus":
                    fseq = comp_seq(samp_fa[qseqid])
                else:
                    fseq = samp_fa[qseqid]
                out_fasta.write(f">{qseqid}|{sseqid}\n{fseq}\n")
                if sseqid not in out_tab_res:
                    out_tab_res[sseqid] = 1.
                else:
                    out_tab_res[sseqid] += 1.
                nb_seq += 1
        for qseqid, c in out_tab_res.items():
            out_tab.write(f"{qseqid} {c}\n")
    

pool = Pool(30)

resutls = pool.map(check_func, (bl_out for bl_out in glob("./blast_res/*.bl")))
