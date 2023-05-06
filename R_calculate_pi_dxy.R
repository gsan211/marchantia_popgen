
library(tidyr)
library(data.table)
library(windowscanr)
library(dplyr)
library(matrixStats)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)

################################################################################################################################################
################################################################################################################################################
#divergence between subspecies based analysis
df = fread("/plas1/george.sandler/marchantia_popgen_newgenome/adhoc_admix/Mar_both_sbsp_filtered_autosomes_GTS.txt")
data_raw = df


cols = ncol(df)

df$missing_rud <- rowSums(df[,5:19] ==".", na.rm=T)
df$missing_poly <- rowSums(df[,20:cols] ==".", na.rm=T)


df2 = df[df$missing_rud < 8 &  df$missing_poly < 5,]

df2$missing <- rowSums(df2[,5:cols] ==".", na.rm=T)
df2$rud_af <- rowSums(df2[,5:19] ==1, na.rm=T)/(rowSums(df2[,5:19] ==0, na.rm=T) + rowSums(df2[,5:19] ==1, na.rm=T))
df2$poly_af <- rowSums(df2[,20:cols] ==1, na.rm=T)/(rowSums(df2[,20:cols] ==0, na.rm=T) + rowSums(df2[,20:cols] ==1, na.rm=T))




#following pi estimates will differ slightly from estimates reported within species since the input data here require coverage in both species resulting in a slightly different subset of sites

#calculate pi in polymorpha
df2$alt <- rowSums(df2[,20:cols] ==1, na.rm=T)
df2$ref   <- rowSums(df2[,20:cols] ==0, na.rm=T)
df2$poly_pi = 1- (choose(df2$ref, 2)+ choose(df2$alt, 2))/(choose(df2$ref + df2$alt, 2))

#calculate pi in ruderalis
df2$alt_rud <- rowSums(df2[,5:19] ==1, na.rm=T)
df2$ref_rud <- rowSums(df2[,5:19] ==0, na.rm=T)
df2$rud_pi = 1- (choose(df2$ref_rud, 2)+ choose(df2$alt_rud, 2))/(choose(df2$ref_rud + df2$alt_rud, 2))


#calculate d_xy
df2$alt_rud <- rowSums(df2[,5:19] ==1, na.rm=T)
df2$ref_rud <- rowSums(df2[,5:19] ==0, na.rm=T)


df2$d_xy =  (   ( df2$alt_rud  *  df2$ref ) + (df2$ref_rud  *  df2$alt) )/((15 - df2$missing_rud)*(8 - df2$missing_poly))


mean(df2$d_xy)
mean(df2$poly_pi,na.rm=T)
mean(df2$rud_pi,na.rm=T)
