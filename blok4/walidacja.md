Walidacja wyników
Kontrole podstawowe (15.1)
Ilość wierszy nagłówków:4
Liczba ariantów:1543
awk -F'\t' '!/^#/ { print NF }' final.vcf | sort -u wynik:103
grep -v '^#' final.vcf | cut -f3 | sort | uniq -d | wc -l wynik:0

Kontrola narzędziem bcftools


SN      0       number of records:      1543
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 1543
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0

Kontrola spójności z genomem(15.3):
genom = DArT REF: 367
genom = DArT ALT: 338
poza {REF,ALT}:  838
spójnych (REF lub ALT): 45.7%

