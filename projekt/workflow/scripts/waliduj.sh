set -euo pipefail
export LC_ALL=C

VCF=$1; GENOM=$2; PROBES=$3; RAW=$4; FILT=$5; BEST=$6
STAT_MD=$7; WAL_MD=$8; MIN_CONS=$9

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
FAIL=0

N_PROBES=$(grep -c '^>' "$PROBES")
N_HITS=$(cut -f1 "$RAW" | sort -u | wc -l)
N_FILT_ROWS=$(wc -l < "$FILT")
N_FILT_PROBES=$(cut -f1 "$FILT" | sort -u | wc -l)
N_BEST=$(wc -l < "$BEST")
N_BEST_UNIQ=$(cut -f1 "$BEST" | sort -u | wc -l)
N_VCF=$(grep -vc '^#' "$VCF")
N_SAMPLES=$(awk -F'\t' '/^#CHROM/ { print NF - 9; exit }' "$VCF")
STRANDS=$(grep -v '^#' "$VCF" | cut -f8 | sort | uniq -c \
          | awk '{ printf "%s: %s  ", $2, $1 }')

N_META=$(grep -c '^##' "$VCF")
COLS=$(awk -F'\t' '!/^#/ { print NF }' "$VCF" | sort -u)
N_COLSET=$(echo "$COLS" | wc -l)
COLS_TXT=$(echo "$COLS" | paste -sd ',' -)
DUPS=$(grep -v '^#' "$VCF" | cut -f3 | sort | uniq -d | wc -l)

GT_TABLE=$(grep -v '^#' "$VCF" | awk -F'\t' '
  { for (i = 10; i <= NF; i++) n[$i]++ }
  END { for (g in n) printf "%s\t%d\n", g, n[g] }' | sort)
GT_BAD=$(echo "$GT_TABLE" | awk -F'\t' '
  $1 != "0/0" && $1 != "0/1" && $1 != "1/1" && $1 != "./." { b += $2 }
  END { print b + 0 }')

S_COLS=OK;  [ "$N_COLSET" -eq 1 ]      || { S_COLS=BŁĄD;  FAIL=1; echo "BŁĄD: różna liczba kolumn: $COLS_TXT" >&2; }
S_DUPS=OK;  [ "$DUPS" -eq 0 ]          || { S_DUPS=BŁĄD;  FAIL=1; echo "BŁĄD: powtórzone ID: $DUPS" >&2; }
S_GT=OK;    [ "$GT_BAD" -eq 0 ]        || { S_GT=BŁĄD;    FAIL=1; echo "BŁĄD: niedozwolone genotypy: $GT_BAD" >&2; }
S_BEST=OK;  [ "$N_BEST" -eq "$N_BEST_UNIQ" ] || { S_BEST=BŁĄD; FAIL=1; echo "BŁĄD: blast_best ma powtórzone sondy" >&2; }

if bcftools view "$VCF" > /dev/null 2> "$TMP/view.err"; then
  S_VIEW=OK
else
  S_VIEW=BŁĄD; FAIL=1; cat "$TMP/view.err" >&2
fi
bcftools stats "$VCF" > "$TMP/stats.txt" 2>/dev/null || true
SN=$(grep '^SN' "$TMP/stats.txt" | cut -f3,4)
TSTV=$(grep '^TSTV' "$TMP/stats.txt" | cut -f5 | head -n 1)

samtools faidx "$GENOM"
grep -v '^#' "$VCF" | awk -F'\t' '{ print $1 ":" $2 "-" $2 }' > "$TMP/regiony.txt"
samtools faidx -r "$TMP/regiony.txt" "$GENOM" | grep -v '^>' > "$TMP/zasady.txt"

if [ "$(wc -l < "$TMP/regiony.txt")" -ne "$(wc -l < "$TMP/zasady.txt")" ]; then
  echo "BŁĄD: liczba pobranych zasad różni się od liczby wariantów" >&2; exit 1
fi

paste <(grep -v '^#' "$VCF" | cut -f1-5,8) "$TMP/zasady.txt" | awk -F'\t' '
  {
    s = $6; sub(/^STRAND=/, "", s)
    g = toupper($7); r = toupper($4); a = toupper($5)
    n[s]++; N++
    if (g == r)      { er++; ok[s]++ }
    else if (g == a) { ea++; ok[s]++ }
    else             { poza++ }
  }
  END {
    printf "ER=%d\nEA=%d\nPOZA=%d\nPCT=%.1f\n", er, ea, poza, (N ? 100*(er+ea)/N : 0)
    for (s in n) printf "STR %s\t%d\t%d\t%.1f\n", s, n[s], ok[s]+0, 100*(ok[s]+0)/n[s]
  }' > "$TMP/spojnosc.txt"

ER=$(grep '^ER=' "$TMP/spojnosc.txt" | cut -d= -f2)
EA=$(grep '^EA=' "$TMP/spojnosc.txt" | cut -d= -f2)
POZA=$(grep '^POZA=' "$TMP/spojnosc.txt" | cut -d= -f2)
PCT=$(grep '^PCT=' "$TMP/spojnosc.txt" | cut -d= -f2)
STR_ROWS=$(grep '^STR ' "$TMP/spojnosc.txt" | sed 's/^STR //' \
           | awk -F'\t' '{ printf "| %s | %s | %s | %s%% |\n", $1, $2, $3, $4 }')

S_CONS=OK
if awk -v p="$PCT" -v m="$MIN_CONS" 'BEGIN { exit !(p < m) }'; then
  S_CONS=BŁĄD; FAIL=1
  echo "BŁĄD: spójność z genomem ${PCT}% < ${MIN_CONS}% (sprawdź pozycję i korektę nici)" >&2
fi


mkdir -p "$(dirname "$STAT_MD")" "$(dirname "$WAL_MD")"

cat > "$STAT_MD" <<EOF
# Statystyki potoku DArT → VCF

| Etap | Liczba |
| :--- | ---: |
| Sondy wejściowe (probes.fa) | $N_PROBES |
| Sondy z trafieniem BLAST | $N_HITS |
| Wiersze trafień po filtrze | $N_FILT_ROWS |
| Sondy po filtrze | $N_FILT_PROBES |
| Markery po wyborze najlepszego | $N_BEST |
| Warianty w final.vcf | $N_VCF |
| Próbki | $N_SAMPLES |

Orientacja nici: $STRANDS
EOF

cat > "$WAL_MD" <<EOF
# Walidacja final.vcf

## Kontrole podstawowe

| Kontrola | Wynik | Status |
| :--- | :--- | :---: |
| Wiersze nagłówkowe (##) | $N_META | |
| Liczba kolumn w wierszach danych (jedna wartość) | $COLS_TXT | $S_COLS |
| Powtórzone ID (oczekiwane 0) | $DUPS | $S_DUPS |
| Niedozwolone genotypy (oczekiwane 0) | $GT_BAD | $S_GT |
| Jeden wiersz na marker w blast_best | $N_BEST / $N_BEST_UNIQ | $S_BEST |
| bcftools view | wczytuje plik | $S_VIEW |

Rozkład genotypów:

$(echo "$GT_TABLE" | awk -F'\t' '{ printf "- %s: %s\n", $1, $2 }')

## bcftools stats

$(echo "$SN" | awk -F'\t' '{ printf "- %s %s\n", $1, $2 }')
- ts/tv: $TSTV (typowo 1,5 do 3)

## Spójność z genomem (zasada ∈ {REF, ALT})

- genom = DArT REF: $ER
- genom = DArT ALT: $EA
- poza {REF, ALT}: $POZA
- spójnych: **${PCT}%** (próg ${MIN_CONS}%): $S_CONS

| Nić | Markery | Zgodne | Odsetek |
| :--- | ---: | ---: | ---: |
$STR_ROWS
EOF

[ "$FAIL" -eq 0 ] || exit 1