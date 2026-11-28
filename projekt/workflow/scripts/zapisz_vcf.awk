BEGIN { FS = OFS = "\t" }

NR == 1 {
    print "##fileformat=VCFv4.2"
    print "##source=DArT-BLAST-pipeline"
    print "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Orientacja dopasowania BLAST\">"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    for (i = 7; i <= NF; i++) printf "\t%s", $i
    printf "\n"
    next
}

{
    printf "%s\t%s\t%s\t%s\t%s\t.\tPASS\tSTRAND=%s\tGT", $1, $2, $3, $4, $5, $6
    for (i = 7; i <= NF; i++) printf "\t%s", $i
    printf "\n"
}