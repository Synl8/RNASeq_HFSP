import gzip


def read_designs(infile="/home/vaitea/project_ori/sequencing_2/data_curie/analysis/sequences/all_seq.fa"):
    seen = {}
    all_design_type = []
    name_id = 0
    for name, seq in read_fasta(infile):
        if len(name.split("_")) > 1:
            name_ = "_".join(name.split("_")[:-1])
        else:
            name_ = name
        if "DCA_SB" in name_:
            name_ = "DCA_SB"
        if name_ not in seen:
            seen[name_] = name_id
            name_id += 1
        all_design_type += [seen[name_]]
    return all_design_type


def comp_seq(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join([comp[s] for s in seq][::-1])


def enum_fasta(infile):
    "for large fasta files, use enumerator instead"
    name = ""
    seq = ""
    for l in open(infile):
        if l.startswith(">"):
            if seq != "":
                yield (name, seq)
                seq = ""
            name = l[1:].strip()
        else:
            seq += l.strip("\n\r")
    if seq != "":
        yield (name, seq)


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


def read_fastq_bin(infile):
    for l_ in gzip.open(infile, "rb"):
        l = str(l_, "utf-8")
        # name
        if l.startswith("@"):
            name = l.strip()
            name_flag = True
        elif name_flag:
            name_flag = False
            seq = l.strip()
            yield seq


def read_fastq_qual(infile):
    fastq_stream = gzip.open(infile)
    while fastq_stream.readline():
        seq = str(fastq_stream.readline(), "utf-8").strip()
        fastq_stream.readline()
        qual = str(fastq_stream.readline(), "utf-8").strip()
        yield seq, [ord(el) for el in qual]


def read_fastq_qual_txt(infile):
    fastq_stream = open(infile)
    while fastq_stream.readline():
        seq = fastq_stream.readline().strip()
        fastq_stream.readline()
        qual = fastq_stream.readline().strip()
        yield seq, [ord(el) for el in qual]


def read_fastq(infile):
    fastq_stream = open(infile)
    while fastq_stream.readline():
        seq = fastq_stream.readline().strip()
        fastq_stream.readline()
        fastq_stream.readline().strip()
        yield seq


from numpy import array
NUC = "ACUTGK-"
NUC_B = [[1., 0., 0., 0., 0.],
         [0., 1., 0., 0., 0.],
         [0., 0., 1., 0., 0.],
         [0., 0., 1., 0., 0.],
         [0., 0., 0., 1., 0.],
         [0., 0., 0., 0., 0.],
         [0., 0., 0., 0., 0.]]


def encode_seq(seq):
    res = []
    for el in seq:
        try:
            res += NUC_B[NUC.index(el)]
        except:
            print(el, seq)
            break
    return res
