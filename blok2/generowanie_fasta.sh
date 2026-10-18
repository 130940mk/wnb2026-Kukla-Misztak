#!\bin\bash 

awk -F ',' 'NR>=6 && (NR-6)%2==0 {print(">"$1 "\n"$2)}' DArT_SNP.csv > probes.fa

#Wyniki kontrolne:
#bioinfo_base.yml  probes.fa  srodowisko.yml  start_Kukla.sh
#(base) [wcss] michal97@r11ch01b02 ~/wnb2026 > grep -c '>' probes.fa
#2008
#(base) [wcss] michal97@r11ch01b02 ~/wnb2026 > awk 'NR%2==0 { print length($0) }' probes.fa | sort -u
#69
