#!/usr/bin/env bash

blast=./blast_bin/bin/
db_blast="./db/pool_0/db/pool0"

for reads in reads/*.fa; do
    name=$(basename $reads .fa)
    $blast/blastn -query $reads -db $db_blast \
                  -outfmt "6 qseqid sseqid pident qcovs qcovus eval sstrand sseq qseq" \
                  -num_threads 30 -max_target_seqs 1 > blast_res/${name}.bl
done
