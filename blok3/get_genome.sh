#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
    AR=$line
    echo "*** Pobieranie: $AR"
    esearch -db biosample -query $AR < /dev/null \
      | elink -target assembly \
      | efetch -format docsum \
      | xtract -pattern DocumentSummary -element FtpPath_GenBank \
      | xargs -n 1 sh -c 'wget "$0"/*_genomic.fna.gz'
done < $1
