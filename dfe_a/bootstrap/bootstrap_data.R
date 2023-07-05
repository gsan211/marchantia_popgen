
#library(tidyr)
library(data.table)
#library(dplyr)

args = commandArgs(trailingOnly = TRUE)



file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/Mar_all_filtered_", args[1] ,".GTS.txt")
print(file_in)


neut = fread(file_in)
#neut2 = neut[neut$V1 == "chr1",]


neut = neut[neut$V1 == "chr1" | neut$V1 == "chr2" | neut$V1 == "chr3" | neut$V1 == "chr4" | neut$V1 == "chr5" | neut$V1 == "chr6" | neut$V1 == "chr7" | neut$V1 == "chr8" ,]


neut$alt <- rowSums(neut[,5:19] ==1, na.rm=T)
neut$ref <- rowSums(neut[,5:19] ==0, na.rm=T)
neut$missing <- rowSums(neut[,5:19] ==".", na.rm=T)


neut2 = neut[neut$missing == 4,]
neut3 = neut2[sample(nrow(neut2), nrow(neut2), replace=TRUE), ]

sfs0 = data.frame(table(neut3$alt))


sfs0$Var1 = as.numeric(as.character(sfs0$Var1))


##############################################################################################################################################

neut = fread("/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_4fold_degen.GTS.txt")
#neut2 = neut[neut$V1 == "chr1",]


neut = neut[neut$V1 == "chr1" | neut$V1 == "chr2" | neut$V1 == "chr3" | neut$V1 == "chr4" | neut$V1 == "chr5" | neut$V1 == "chr6" | neut$V1 == "chr7" | neut$V1 == "chr8" ,]


neut$alt <- rowSums(neut[,5:19] ==1, na.rm=T)
neut$ref <- rowSums(neut[,5:19] ==0, na.rm=T)
neut$missing <- rowSums(neut[,5:19] ==".", na.rm=T)


neut2 = neut[neut$missing == 4,]
neut3 = neut2[sample(nrow(neut2), nrow(neut2), replace=TRUE), ]


sfs4 = data.frame(table(neut3$alt))


sfs4$Var1 = as.numeric(as.character(sfs4$Var1))




##############################################################################################################################################
##############################################################################################################################################
Var1 = c(0:11)
df = data.frame(Var1) #create vector of appropriate length to check for missing SFS bins


sfs_unfolded = merge(sfs0, sfs4, by = "Var1", all=T) #check for missing categories here!!
sfs_unfolded2 = merge(sfs_unfolded, df, by = "Var1", all=T) #merge with dataframe of correct length, creates NA's for any missing SFS bins, replace them with 0
sfs_unfolded2[is.na(sfs_unfolded2)] <- 0
sfs_unfolded2$Var1 = NULL

sfs_unfolded3 = t(sfs_unfolded2)


write.table(sfs_unfolded3, "/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/bootstrapped_sfs", sep = "\t", col.names =F, quote = F, row.names =F)

##############################################################################################################################################
##############################################################################################################################################
 
div_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/divergence_file_", args[1] ,".txt")
print(div_in)

div = fread(div_in)
#div = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/divergence_file_3prime_UTR.txt")

#use binomial sampling to generate bootstrapped divergence matrix

#first calculate "success" rate for each trial based on known data
prob_1 = as.numeric(div[1,3]/div[1,2])
prob_0 = as.numeric(div[2,3]/div[2,2])

#then count up "#successes" based on data
n_div1 = sum(rbinom(div[[1,2]],1,prob_1))
n_div0 = sum(rbinom(div[[2,2]],1,prob_0))

#copy existing df just easy to preserve dfe input structure and sub in the new smapled estimates
div2 = div
div2[1,3] = n_div1
div2[2,3] = n_div0

write.table(div2, "/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/bootstrap/bootstrapped_divergence", sep = "\t", col.names =F, quote = F, row.names =F)














