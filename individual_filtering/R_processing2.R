#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

print(args[1])

library(tidyr)
library(data.table)
library(windowscanr)
library(dplyr)
library(matrixStats)
library(GenomicRanges)



infile = paste("/plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/windows_analysis_output/",args[1],"_allwindows.txt",sep="")

#df = fread("/plas1/george.sandler/marchantia_popgen/filtering/polymorpha_individual/windows_analysis_output/Mar_E_allwindows.txt", header=T)


df = fread(infile, header =F)

df$qual = df$V4/df$V5



df2 = df[df$qual > 0.00572,] #choose bad windows with ~2 or more bad SNPs 


dt <- makeGRangesFromDataFrame(df2, seqnames.field = "V1", start.field = "V2", end.field = "V3")
r <- reduce(dt) #merge overlapping ranges
dc <- data.frame(r) #convert back to df
dc$width <- NULL
dc$strand <- NULL


outfile = paste("/plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/windows_analysis_output/",args[1],"_bad_windows.txt",sep="")
write.table(dc, outfile,row.names = F, col.names = F, quote = F, sep = "\t")
