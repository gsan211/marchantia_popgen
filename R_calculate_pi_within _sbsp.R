
library(tidyr)
library(data.table)
library(windowscanr)
library(dplyr)
library(matrixStats)
library(GenomicRanges)

df = fread("/plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround.GTS.txt")

#header order:    Mar_beta Mar_D   Mar_E   Mar_G   Mar_H   Mar_J   Mar_K   Mar_L   Mar_M   Mar_N   Mar_O   Mar_S   Mar_T   Mar_U   Mar_Y

#ellice   Mar_J   Mar_K   Mar_L  Mar_N   Mar_O   Mar_S   Mar_T   Mar_U
#df = df[,c(1:4,10,11,12,14,15,16,17,18)]        #select ellice only
tot = ncol(df)
df$alt <- rowSums(df[,5:tot] ==1, na.rm=T)
df$ref <- rowSums(df[,5:tot] ==0, na.rm=T)
df$missing <- rowSums(df[,5:tot] ==".", na.rm=T)


#calculate denominator and numerator of pi for each site
df$deno_pi =  choose(df$ref + df$alt, 2)
df$nume_pi = choose(df$ref, 2)+ choose(df$alt, 2)


#select major assembled chromosomes
df = df[df$V1 == "chr1" | df$V1 == "chr2" | df$V1 == "chr3" | df$V1 == "chr4" | df$V1 == "chr5" | df$V1 == "chr6" | df$V1 == "chr7" | df$V1 == "chr8" ,]

df$pi = 1- df$nume_pi/df$deno_pi



#calculate genome-wide average pi by summing pi across all variant sites
sum(df$pi, na.rm=T) /nrow(df)



######################################################################################################################################
######################################################################################################################################

#pairwise differences between all pairs of samples
pair = data.frame(combn(5:19, 2))
mat = NULL
for(i in 1:105){

  a = pair[1,i]
  b = pair[2,i]
  vec = c(1:4,a,b)
  ab <-   df %>% select(1:5,a,b)
  #ab = df[,1:7]
  ab$alt <- rowSums(ab[,5:6] ==1, na.rm=T)
  ab$ref <- rowSums(ab[,5:6] ==0, na.rm=T)
  #ab$len = ab$V4 - ab$V2
  ab = ab[(ab$ref + ab$alt) > 1,]
  ab$deno_pi =  choose(ab$ref + ab$alt, 2)
  ab$nume_pi = choose(ab$ref, 2)+ choose(ab$alt, 2)

  ab$pi = 1- ab$nume_pi/ab$deno_pi

  
  pi = mean(ab$pi)

  line = data.frame(noquote(paste(a,b,pi)))
  mat = rbind(line,mat)
}

write.table(mat, "/plas1/george.sandler/marchantia_popgen_newgenome/pairwise_pi.txt", col.names = F, row.names =F, quote = F)
