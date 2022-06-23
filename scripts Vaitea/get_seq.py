"""Read the fastq for all reads and create a fasta file
"""

from src.seq_io import read_fastq
from glob import glob
import re

for in_fastq in glob("./fastq/*.notCombined_1.fastq"):
    base_name = re.match("./fastq/(.+)\.notCombined_1\.fastq", in_fastq).groups()[0]
    with open(f"reads/{base_name}_1.fa", "w") as out:
        for si, seq in enumerate(read_fastq(in_fastq)):
            out.write(f">{si}\n{seq}\n")
