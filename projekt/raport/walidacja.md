# Walidacja final.vcf

## 1. Kontrole podstawowe (15.1)

| Kontrola | Polecenie | Wynik |
|---|---|---|
| wiersze nagłówkowe | `grep -c '^##' final.vcf` | 4 |
| warianty | `grep -vc '^#' final.vcf` | 1543 |
| liczba kolumn | `awk -F'\t' '!/^#/ {print NF}' final.vcf \| sort -u` | 103 |
| duplikaty ID | `grep -v '^#' final.vcf \| cut -f3 \| sort \| uniq -d \| wc -l` | 0 |

## 2. bcftools (15.2)

```
SN  0  number of samples:  94
SN  0  number of records:  1543
SN  0  number of SNPs:     1543
SN  0  number of indels:   0
```

## 3. Spójność z genomem (15.3)

| Kategoria | Liczba |
|---|---|
| genom = DArT REF | 367 |
| genom = DArT ALT | 338 |
| poza {REF, ALT} | 838 |
| spójnych | 45,7% |
