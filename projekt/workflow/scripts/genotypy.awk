BEGIN { FS = ","; OFS = "\t" }

NR == header_row {
    nsamp = NF - sample_col + 1
    printf "ID"
    for (i = sample_col; i <= NF; i++) printf OFS $i
    printf "\n"
    next
}

NR < data_start { next }

{
    start = NF - nsamp + 1

    if (prev == "") {
        prev = $1
        for (i = start; i <= NF; i++) ref[i - start] = $i
        next
    }

    printf "%s", $1
    for (i = start; i <= NF; i++) {
        r = ref[i - start]
        a = $i
        if (r == "-" || a == "-") gt = "./."
        else if (r == 1 && a == 0) gt = "0/0"
        else if (r == 0 && a == 1) gt = "1/1"
        else if (r == 1 && a == 1) gt = "0/1"
        else gt = "./."
        printf OFS gt
    }
    printf "\n"
    prev = ""
}