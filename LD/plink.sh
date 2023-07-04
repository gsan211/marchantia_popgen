cd "/plas1/george.sandler/marchantia_popgen/plink/"




########################## randomly sample sites
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 
do

grep -v "#" "/plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround_autosomes_variant_noaltfixed.vcf" | grep -e ${chr} | shuf -n 5000 > /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom1.vcf


grep -e "#" "/plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround_autosomes_variant_noaltfixed.vcf" > "/plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_header.vcf"

cat "/plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_header.vcf" /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom1.vcf >  /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom.vcf

rm  /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom1.vcf
rm "/plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_header.vcf"

done


########################## 


#LD calculations
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 
do

"/plas1/george.sandler/apps/plink" --make-bed  --vcf /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom.vcf --out /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --double-id  --allow-extra-chr 


# --maf 0.001
"/plas1/george.sandler/apps/plink" --bfile /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --recode --tab  --out /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --allow-extra-chr


#"/plas1/george.sandler/apps/plink" --made bed --file /plas1/george.sandler/marchantia_popgen/plink/plink --maf 0.01 -out /plas1/george.sandler/marchantia_popgen/plink/plink2

"/plas1/george.sandler/apps/plink" --file /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --r  --out /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom --allow-extra-chr  --with-freqs --inter-chr 

rm /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr?_5krandom.log
rm /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr?_5krandom.nosex
done


#

#"/plas1/george.sandler/apps/plink" --file /plas1/george.sandler/marchantia_popgen/plink/plink  --freq --out /plas1/george.sandler/marchantia_popgen/plink/UV_allind --allow-extra-chr --inter-chr

#####################


#merge plink LD
cat /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr?_5krandom.ld | grep -v "nan" > /plas1/george.sandler/marchantia_popgen_newgenome/ld/Mar_all_filtered_all_5krandom.ld

(grep -e "CHR" /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr1_5krandom.ld ; sed  's/.*CHR.*//g' /plas1/george.sandler/marchantia_popgen_newgenome/ld/Mar_all_filtered_all_5krandom.ld) | grep . > /plas1/george.sandler/marchantia_popgen_newgenome/ld/Mar_all_filtered_all_5krandom_final.ld

rm "/plas1/george.sandler/marchantia_popgen_newgenome/ld/Mar_all_filtered_all_5krandom.ld"

#################################

#sex chromosomes


#LD calculations
for chr in U V
do

"/plas1/george.sandler/apps/plink" --make-bed  --vcf /plas1/george.sandler/marchantia_popgen_newgenome/calls/UV/Mar_all_${chr}_secondfiltered_variant_noaltfixed.vcf --out /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --double-id  --allow-extra-chr 


# --maf 0.001
"/plas1/george.sandler/apps/plink" --bfile /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --recode --tab  --out /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --allow-extra-chr


#"/plas1/george.sandler/apps/plink" --made bed --file /plas1/george.sandler/marchantia_popgen/plink/plink --maf 0.01 -out /plas1/george.sandler/marchantia_popgen/plink/plink2

"/plas1/george.sandler/apps/plink" --file /plas1/george.sandler/marchantia_popgen_newgenome/plink/plink --r  --out /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_${chr}_5krandom --allow-extra-chr  --with-freqs --inter-chr 

rm /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr_?_5krandom.log
rm /plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_chr_?_5krandom.nosex
done












##########LD pruning
#first need to add site info to annotation field
#bgzip
#tabix

/plas1/george.sandler/apps/bcftools-1.10.2/bin/bcftools annotate --set-id +'%CHROM\-%POS' /plas1/george.sandler/marchantia_popgen/structure/Mar_all_perfectSNPs_shared.vcf > "/plas1/george.sandler/marchantia_popgen/structure/Mar_all_perfectSNPs_shared_siteids.vcf"

cd "/plas1/george.sandler/marchantia_popgen/plink/"
#prep vcf to plink
"/plas1/george.sandler/apps/plink" --make-bed  --vcf "/plas1/george.sandler/marchantia_popgen/structure/Mar_all_perfectSNPs_shared_siteids.vcf" --out /plas1/george.sandler/marchantia_popgen/plink/Mar_all_shared_perfectSNPs --double-id  --allow-extra-chr

"/plas1/george.sandler/apps/plink" --bfile /plas1/george.sandler/marchantia_popgen/plink/Mar_all_shared_perfectSNPs --recode --tab  --out /plas1/george.sandler/marchantia_popgen/plink/Mar_all_shared_perfectSNPs --allow-extra-chr

#pruning
cd /plas1/george.sandler/marchantia_popgen/structure/
"/plas1/george.sandler/apps/plink" --file /plas1/george.sandler/marchantia_popgen/plink/Mar_all_shared_perfectSNPs --indep 100 5 1.5 --allow-extra-chr


awk -F '-' 'BEGIN {OFS="\t"}; {print $1,$2}' "/plas1/george.sandler/marchantia_popgen/structure/plink.prune.in" > "/plas1/george.sandler/marchantia_popgen/structure/LD_pruned_perfect_SNPs_shared.list"



#done

R
library(data.table)

df = fread("/plas1/george.sandler/marchantia_popgen_newgenome/plink/Mar_all_filtered_secondround_autosomes_4fold_quart4.ld")
df = na.omit(df)
df = df[df$CHR_A == df$CHR_B,]

df$dist = df[,6]-df[,2]
df$R2 = df$R*df$R
mod = smooth.spline(df$dist,df$R2)
plot(mod)

#plot(mod, xlim = c(0,1000000))
plot(mod, xlim = c(0,25000))



df2 = df[df$CHR_A != "Chr_Y_A" & df$CHR_A != "Chr_Y_B" ,]
df2 = df2[df2$CHR_A != "scaffold_18" & df2$CHR_A != "scaffold_17" ,]

mod2 = smooth.spline(df2$dist,df2$R2)
plot(mod2, xlim = c(0,25000))


df3 = df[df$CHR_A == "Chr_Y_A" | df$CHR_A == "Chr_Y_B" |df$CHR_A == "scaffold_18" | df$CHR_A == "scaffold_17"  ,]

mod3 = smooth.spline(df3$dist,df3$R2)
plot(mod3, xlim = c(0,25000))
lines(mod3, col='red', lwd=2)



df$bin = (cut(df$dist, c(0,10,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,5000,10000,25000,100000,Inf)))
df2 = aggregate(df$R2~ df$bin, df, mean)

barplot(df2[,2],col="black", main="Sum of D (netLD) LOF Capsella max5" )
####################
#interchromosomal comparisons


df = fread("/plas1/george.sandler/marchantia_popgen/plink/Mar_sbsp_polymorpha_perfectSNPs_3.5Krandom.ld")
df = na.omit(df)
df = df[df$CHR_A == df$CHR_B,]
df$R2 = df$R*df$R


df$inter <- ifelse(df$CHR_A == df$CHR_B, 0, 1)
