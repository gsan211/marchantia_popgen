bcftoolsg=/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools
samtoolsg=/plas1/george.sandler/apps/samtools-1.8/samtools


#cd "/plas1/george.sandler/apps/snpEff/"
#java -jar snpEff.jar build -gff3 -v marchantia

cd "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/"
#prepare VCF for full annotation:
grep -v "#" "/plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround.vcf" | awk -F"\t" '{$5="N"; print}' OFS="\t"   > /plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround_edited1.vcf


grep -e "#" "/plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround.vcf" > "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/tmpHeader.vcf"

cat "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/tmpHeader.vcf" /plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround_edited1.vcf > /plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround_edited.vcf

rm "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/tmpHeader.vcf"
rm /plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround_edited1.vcf

############### run snpeff

java -jar "/plas1/george.sandler/apps/snpEff/snpEff.jar" -c /plas1/george.sandler/apps/snpEff/snpEff.config  -v marchantia_new /plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround.vcf    > /plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround.ann.vcf

############# grep out 0-fold and 4-fold sites

grep -e "#" -e "synonymous" "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround_edited.ann.vcf" | grep -Ev 'missense|HIGH|ERROR'  > "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_4fold_degen.ann.vcf" & 

grep -e "#" -e "missense" "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround_edited.ann.vcf" | grep -Ev "synonymous|HIGH|ERROR" > "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_0fold_degen.ann.vcf" &

#count # of matches per line, to check that I'm actually picking sites with the right degeneracy
#grep -o -n 'synonymous' /plas1/george.sandler/marchantia_popgen/snpeff/Mar_polymorpha_filtered_secondround_4fold_degen.ann.vcf2 | cut -d : -f 1 | uniq -c | less -S 

################################################################################################################################################
################################################################################################################################################
############# grep out other site types

#5_prime_UTR
grep -e "#" -e "5_prime" "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround_edited.ann.vcf" | grep -Ev 'missense|HIGH|ERROR|synonymous|intron|3_prime|splice'  > "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_5primeUTR.ann.vcf" &

#3_prime_UTR
grep -e "#" -e "3_prime" "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround_edited.ann.vcf" | grep -Ev 'missense|HIGH|ERROR|synonymous|intron|5_prime|splice'  > "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_3primeUTR.ann.vcf" &


#splice site intron
grep -e "#" -e "splice" "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround_edited.ann.vcf" | grep -Ev 'missense|LOW|MODERATE|ERROR|synonymous|3_prime|5_prime'  > "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_all_splice_site_intron.ann.vcf" &

#all intron
grep -e "#" -e "intron" "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround_edited.ann.vcf" | grep -Ev 'missense|HIGH|ERROR|synonymous|3_prime|5_prime'  > "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_allintron.ann.vcf" &

#all_noncoding
grep -e "#" -e "intergenic" "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_secondround_edited.ann.vcf" | grep -Ev 'missense|HIGH|ERROR|synonymous|3_prime|5_prime_intron|splice'  > "/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_allintergenic.ann.vcf" &


##############
#query VCF to more usable text file. Can then calculate diversity stats with other R scripts etc. 

for ind in Mar_all_filtered_5primeUTR Mar_all_filtered_3primeUTR Mar_all_filtered_all_splice_site_intron Mar_all_filtered_allintron Mar_all_filtered_allintergenic
do

bcftools query /plas1/george.sandler/marchantia_popgen_newgenome/snpeff/${ind}.ann.vcf -f "%CHROM\t%POS\t%ALT\t[\t%GT]\n" -o /plas1/george.sandler/marchantia_popgen_newgenome/snpeff/${ind}.GTS.txt


done
