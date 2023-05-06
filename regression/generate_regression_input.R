
library(tidyr)
library(data.table)
library(windowscanr)
library(dplyr)
library(matrixStats)
library(GenomicRanges)
library(glmnet)  


       
########################################################################################################################################
########################################################################################################################################



cod = fread("/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_CDS_intersect.txt")
cod[,3:6] = NULL
cod$id = paste(cod$V1, cod$V2)

cod$V1 = NULL
cod$V2 = NULL
#need to aggregate by window start to merge counts of seperate genes from bedtools output file

cod2 = aggregate(.~id, cod, FUN=sum)
cod2$frac_cds = cod2$V7/10000    #dont forget to edit based on window size !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cod2$V7 = NULL  

gc = fread("/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_GC_content.txt")   
gc[,c(3,4,6:12)] = NULL
gc$id = paste(gc$'#1_usercol',gc$'2_usercol' )
gc$'#1_usercol' = NULL
gc$'2_usercol' = NULL



codgc = merge(cod2, gc, by = "id", all = T)


########################################################################################################################################
utr = fread("/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_UTR_intersect.txt")
utr[,3:6] = NULL
utr$id = paste(utr$V1, utr$V2)

utr$V1 = NULL
utr$V2 = NULL
#need to aggregate by window start to merge counts of seperate genes from bedtools output file

utr2 = aggregate(.~id, utr, FUN=sum)
utr2$frac_utr = utr2$V7/10000    #dont forget to edit based on window size !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
utr2$V7 = NULL  




cod3 = merge(utr2, codgc, by = "id", all = T)


########################################################################################################################################
########################################################################################################################################
entr = fread("/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_entropy_nomissing_edited.txt")
entr$V3 = NULL
entr$V4 = NULL
entr$V6 = NULL

entr$id = paste(entr$V1, entr$V2)
entr$V1 = NULL
entr$V2 = NULL

entr$entropy = entr$V5
entr$V5 = NULL

ling = fread("/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_lingcomplex_nomissing_edited.txt")
ling$V3 = NULL
ling$V4 = NULL
ling$V6 = NULL

ling$id = paste(ling$V1, ling$V2)
ling$V1 = NULL
ling$V2 = NULL

ling$ling_complex = ling$V5
ling$V5 = NULL

nessie = merge(entr, ling, by ="id", all = T)


############################################################################

#read in pi calculated in 10kb windows


win = fread( "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10kb_win_output_initial_autosomes.txt", header = T)


win2 = win[win$pi_n > 1000,]
win2$winpi = win2$pi_summ/win2$pi_n
win2$winpi4 = win2$degenpi_summ/win2$degenpi_n

win2$id = paste(win2$V1, win2$win_start)

win2$win_start = NULL
win2$win_end= NULL
win2$degenpi_n= NULL
win2$degenpi_summ= NULL
win2$winpi4= NULL
win2$pi_summ = NULL

win3 = merge(win2, cod3, by = "id", all=T)


#win3[["frac"]][is.na(win3[["frac"]])] <- 0 #not necessary with bedtools intersect wao as 0's counted



############################################################################
win4 = merge(win3, nessie, by = "id")


########################################################################################################################################
########################################################################################################################################
quart1 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_quart1_intersect.txt")
quart1[,3:6] = NULL
quart1$id = paste(quart1$V1, quart1$V2)
quart1$V1 = NULL
quart1$V2 = NULL
quart1_2 = aggregate(.~id, quart1, FUN=sum)
quart1= quart1_2[quart1_2$V7 > 1500,]

winq1 = win4[win4$id %in% quart1$id,]
winq1$exp = 1


quart2 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_quart2_intersect.txt")
quart2[,3:6] = NULL
quart2$id = paste(quart2$V1, quart2$V2)
quart2$V1 = NULL
quart2$V2 = NULL
quart2_2 = aggregate(.~id, quart2, FUN=sum)
quart2 = quart2_2[quart2_2$V7 > 1500,]

winq2 = win4[win4$id %in% quart2$id,]
winq2$exp = 2


quart3 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_quart3_intersect.txt")
quart3[,3:6] = NULL
quart3$id = paste(quart3$V1, quart3$V2)
quart3$V1 = NULL
quart3$V2 = NULL
quart3_2 = aggregate(.~id, quart3, FUN=sum)
quart3 = quart3_2[quart3_2$V7 > 1500,]

winq3 = win4[win4$id %in% quart3$id,]
winq3$exp = 3


quart4 = fread("/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/10k_windows_quart4_intersect.txt")
quart4[,3:6] = NULL
quart4$id = paste(quart4$V1, quart4$V2)
quart4$V1 = NULL
quart4$V2 = NULL
quart4_2 = aggregate(.~id, quart4, FUN=sum)
quart4 = quart4_2[quart4_2$V7 > 1500,]

winq4 = win4[win4$id %in% quart4$id,]
winq4$exp = 4


###########
#genes with no expression, non-coding regions and windows where less than 1.5k bases overlap gene
quart0 = rbind(quart1, quart2, quart3, quart4)


winq0 = subset(win4, !(id %in% quart0$id))
winq0$exp = 0



win5 = rbind(winq0 ,winq1, winq2, winq3, winq4)



#remove windows that overlap more than one gene with expression
win5 = win5[!duplicated(win5$id) & !duplicated(win5$id, fromLast = TRUE),]   


############################################################################
#divergence 

ddiv = fread("/plas1/george.sandler/marchantia_popgen_newgenome/Mar_table_poly_d_xy_10kb_allchrom.txt")

ddiv$id = paste(ddiv$V1, ddiv$win_start)
ddiv[,c(1:8)] = NULL


div2 = ddiv[ddiv$d_xy_n > 100,]

win6 = merge(win5, div2, by = "id", all=T)
win6$frac_GC = win6[,8] 
win6[,8]  = NULL

write.table(win6, "/plas1/george.sandler/marchantia_popgen_newgenome/bedtools/regression_input.txt", col.names = T, row.names =F, quote = F)
