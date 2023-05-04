#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

print(args[1])

library(tidyr)
library(data.table)
library(windowscanr)
library(dplyr)
library(matrixStats)

infile = paste("/plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/Mar_bcftools_AD_",args[1],".txt",sep="")
#df = fread("/plas1/george.sandler/marchantia_popgen/filtering/polymorpha_individual/Mar_sbsp_polymorpha_bcftools_AD_Mar_E_old.txt")

df = fread(infile, header =F)

df = df[df$V5 !=".",]


df2= separate(df,V5, c("ref","alt1", "alt2","alt3"),sep = ",")

df2$ref = as.numeric(df2$ref)
df2$alt1 = as.numeric(df2$alt1)
df2$alt2 = as.numeric(df2$alt2)
df2$alt3 = as.numeric(df2$alt3)

df2$alt1[is.na(df2$alt1)] <- 0
df2$alt2[is.na(df2$alt2)] <- 0
df2$alt3[is.na(df2$alt3)] <- 0


df2$alt = df2$alt1+df2$alt2+df2$alt3
df2$ratio = df2$alt/(df2$alt+df2$ref)
df2$bad <- ifelse(df2$ratio > 0.95 | df2$ratio < 0.05, 0, 1)
df2$bad[is.na(df2$bad)] <- 0



outfile = paste("/plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/pileup_sep/",args[1],"_sep.txt",sep="")
write.table(df2, outfile,row.names = F, col.names = T, quote = F, sep = "\t")