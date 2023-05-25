library(data.table)
library(tidyr)
library(dplyr)


df = fread("/plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_unfiltered_variant.gts.txt")

#structure of file is chromosome, site, sample genotypes (0/1/.)
setnames(df, old = c('V1','V2','V4', 'V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15','V16','V17','V18'), new = c('chrom','pos','Mar_beta','Mar_D','Mar_E','Mar_G','Mar_H','Mar_J','Mar_K','Mar_L','Mar_M','Mar_N','Mar_O','Mar_S','Mar_T','Mar_U','Mar_Y'))



#to match hmmIBD requirements
df2 = apply(df[,4:18], 2, function(y) as.numeric(gsub("\\.", "-1", y)))
dfchr = apply(df[,1:2], 2, function(y) as.numeric(gsub("chr", "", y)))

df3 = cbind(dfchr,df2)

df4 = na.omit(df3)


write.table(df4, "/plas1/george.sandler/marchantia_popgen_newgenome/IBD/Mar_all_unfiltered_variant_autosomes_hmmibd_input.txt",sep = "\t", col.names =T, row.names = F, quote = F)

