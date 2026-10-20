blastn \
    -task blastn \
    -db db/genome \
    -query probes.fa \
    -out blast_raw.tsv \
    -evalue 1e-10 \
    -max_target_seqs 1 \
    -num_threads 4 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
