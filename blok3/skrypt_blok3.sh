#!/bin/bash
set -e

# 1) UID złożenia (w środowisku biopython)
chmod +x find_biosample.py
./find_biosample.py
./find_biosample.py ASSEMBLY_UID

# 2) pobranie genomu (w środowisku edirect)
chmod +x get_genome.sh
./get_genome.sh bio_ids.txt

# 3) przygotowanie pliku genomu
gunzip -c GCA_*_genomic.fna.gz > genome.fa
mkdir -p ~/projekt-dart/data/ref
mv genome.fa ~/projekt-dart/data/ref/

# 4) ekstrakcja sond do formatu fasta
awk -F',' 'NR>=6 && (NR-6)%2==0 { print ">"$1"\n"$2 }' DArT_SNP.csv > probes.fa

# 5) ekstracja informacji o SNP
awk -F',' -v OFS='\t' 'NR>=7 && (NR-7)%2==0 { print $1, $11 }' DArT_SNP.csv > dart_info_oneline.tsv

# 6) budowa bazy (w trybie interaktywnym i środowisku blast)
mkdir -p db
makeblastdb -in data/ref/genome.fa -dbtype nucl -out db/genome

# 7) mapowanie sond
chmod +x mapowanie.sh
./mapowanie.sh

# 8) filtrowanie trafień i wybór najlepszego dpoasowania
awk -F'\t' -v OFS='\t' '$3 >= 96 && $4 >= 68' blast_raw.tsv > blast_filtered.tsv
sort -k1,1 -k12,12gr blast_filtered.tsv \
  | awk -F'\t' '!widziane[$1]++' > blast_best.tsv

# 9) wyznaczenie orientacji nici
awk -F'\t' -v OFS='\t' '{ print $0, ($9 < $10 ? "+" : "-") }' blast_best.tsv > dart_blast_filtered.tsv

