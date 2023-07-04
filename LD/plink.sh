cd "/plas1/george.sandler/marchantia_popgen/plink/"




########################## randomly sample sites
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 
do

#randomly sample 5k sites not from header, merge with header and prepare new vcf for input to plink, loop over all main autosomes seprately
grep -v "#" "/plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround_autosomes_variant_noaltfixed.vcf" | grep -e ${chr} | shuf -n 5000 > /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom1.vcf


grep -e "#" "/plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround_autosomes_variant_noaltfixed.vcf" > "/plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_header.vcf"

cat "/plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_header.vcf" /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom1.vcf >  /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom.vcf

rm  /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom1.vcf
rm "/plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_header.vcf"

done


########################## 

#run plink on randomly sampled vcf's created in previuos step
#LD calculations
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 
do

"/plas1/george.sandler/apps/plink" --make-bed  --vcf /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom.vcf --out /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --double-id  --allow-extra-chr 

"/plas1/george.sandler/apps/plink" --bfile /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --recode --tab  --out /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --allow-extra-chr

"/plas1/george.sandler/apps/plink" --file /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --r  --out /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom --allow-extra-chr  --with-freqs --inter-chr 

rm /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr?_5krandom.log
rm /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr?_5krandom.nosex
done


#####################


#merge plink LD output to singlefile
cat /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr?_5krandom.ld | grep -v "nan" > /plas1/george.sandler/marchantia_popgen_newgenome/ld/Mar_all_filtered_all_5krandom.ld

(grep -e "CHR" /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr1_5krandom.ld ; sed  's/.*CHR.*//g' /plas1/george.sandler/marchantia_popgen_newgenome/ld/Mar_all_filtered_all_5krandom.ld) | grep . > /plas1/george.sandler/marchantia_popgen_newgenome/ld/Mar_all_filtered_all_5krandom_final.ld

rm "/plas1/george.sandler/marchantia_popgen_newgenome/ld/Mar_all_filtered_all_5krandom.ld"

#################################

#sex chromosomes


#LD calculations
for chr in U V
do

"/plas1/george.sandler/apps/plink" --make-bed  --vcf /plas1/george.sandler/marchantia_popgen_newgenome/calls/UV/Mar_all_${chr}_secondfiltered_variant_noaltfixed.vcf --out /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --double-id  --allow-extra-chr 


"/plas1/george.sandler/apps/plink" --file /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --r  --out /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom --allow-extra-chr  --with-freqs --inter-chr 

rm /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr_?_5krandom.log
rm /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr_?_5krandom.nosex
done


