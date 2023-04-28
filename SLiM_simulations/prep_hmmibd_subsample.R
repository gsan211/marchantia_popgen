library(data.table)
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)




print(args[1])

file_in = paste0(args[1],".vcf")
file_out = paste0("hmmibd/",args[1],"_input.txt")

print(file_in)

print(file_out)
#file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/admix_subpop/slim_subpop_ad65k_10x_",i,".vcf")
#file_out = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/admix_subpop/hmmibd/slim_subpop_ad65k_10x_input",i,".txt")

df = fread(file_in)


df2 = apply(df[,10:ncol(df)], 2, function(y) as.numeric(gsub("\\|.", "", y)))
dfchr = apply(df[,1:2], 2, function(y) as.numeric(gsub("chr", "", y)))
colnames(dfchr) <- c("chrom","pos")


df3 = cbind(dfchr,df2)
df3 = data.frame(df3)


df3 = df3[sample(nrow(df3), nrow(df2)/1.1 ,replace = FALSE), ]

df3 <- df3[order(df3$pos),]

write.table(df3, file_out,sep = "\t", col.names =T, row.names = F, quote = F)



