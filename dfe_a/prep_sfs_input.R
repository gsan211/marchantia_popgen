
library(tidyr)
library(data.table)
library(dplyr)

#read in file from specific genomic feature here, generated from SNPeff + BCFtools
seln = fread("/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/Mar_all_filtered_all_splice_site_intron.GTS.txt")


seln = seln[seln$V1 == "chr1" | seln$V1 == "chr2" | seln$V1 == "chr3" | seln$V1 == "chr4" | seln$V1 == "chr5" | seln$V1 == "chr6" | seln$V1 == "chr7" | seln$V1 == "chr8" ,]


seln$alt <- rowSums(seln[,5:19] ==1, na.rm=T)
seln$ref <- rowSums(seln[,5:19] ==0, na.rm=T)
seln$missing <- rowSums(seln[,5:19] ==".", na.rm=T)


seln2 = seln[seln$missing == 4,]


sfs0 = data.frame(table(seln2$alt))
sfs0$Var1 = as.numeric(as.character(sfs0$Var1))


##############################################################################################################################################

neut = fread("/plas1/george.sandler/marchantia_popgen_newgenome/snpeff/Mar_all_filtered_4fold_degen.GTS.txt")



neut = neut[neut$V1 == "chr1" | neut$V1 == "chr2" | neut$V1 == "chr3" | neut$V1 == "chr4" | neut$V1 == "chr5" | neut$V1 == "chr6" | neut$V1 == "chr7" | neut$V1 == "chr8" ,]


neut$alt <- rowSums(neut[,5:19] ==1, na.rm=T)
neut$ref <- rowSums(neut[,5:19] ==0, na.rm=T)
neut$missing <- rowSums(neut[,5:19] ==".", na.rm=T)


neut2 = neut[neut$missing == 4,]


sfs4 = data.frame(table(neut2$alt))
sfs4$Var1 = as.numeric(as.character(sfs4$Var1))

##############################################################################################################################################

sfs_unfolded = merge(sfs0, sfs4, by = "Var1", all=T)
sfs_unfolded$Var1 = NULL

sfs_unfolded2 = t(sfs_unfolded)


write.table(sfs_unfolded2, "/plas1/george.sandler/marchantia_popgen_newgenome/dfe_a/expanded_categories/SFS_unfolded_input_splice_site_intron.txt", sep = "\t", col.names =F, quote = F, row.names =F)
