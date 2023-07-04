library(data.table)
library(ggplot2)



df = fread("/plas1/george.sandler/marchantia_popgen_newgenome/ld/Mar_all_filtered_all_5krandom_final.ld")
df = na.omit(df)
df = df[df$CHR_A == df$CHR_B,]

df$dist = df[,6]-df[,2]
df$dist2 = df$dist/1000
df$R2 = df$R*df$R
ld = df

#calculate pi
df = fread("/plas1/george.sandler/marchantia_popgen_newgenome/calls/Mar_all_filtered_secondround.GTS.txt")

tot = ncol(df)
df$alt <- rowSums(df[,5:tot] ==1, na.rm=T)
df$ref <- rowSums(df[,5:tot] ==0, na.rm=T)
df$missing <- rowSums(df[,5:tot] ==".", na.rm=T)
df$deno_pi =  choose(df$ref + df$alt, 2)
df$nume_pi = choose(df$ref, 2)+ choose(df$alt, 2)

df = df[df$V1 == "chr1" | df$V1 == "chr2" | df$V1 == "chr3" | df$V1 == "chr4" | df$V1 == "chr5" | df$V1 == "chr6" | df$V1 == "chr7" | df$V1 == "chr8" ,]

df$pi = 1- df$nume_pi/df$deno_pi


#calculate genome-wide average pi by summing pi across all variant sites, and divide by the number of variant sites + length of blocks of invariant sites
sum(df$pi, na.rm=T) /nrow(df)



df2 = aggregate(df$pi,by=list(df$V1), FUN=mean)
colnames(df2) <- c("chr", "pi")


ld1 = ld[ld$dist < 1000,]
ld1a = aggregate(ld1$R2,by=list(ld1$CHR_A), FUN=mean)
ld1a$dist = "<1000"


ld2 = ld[ld$dist > 1000 & ld$dist < 10000,]
ld2a = aggregate(ld2$R2,by=list(ld2$CHR_A), FUN=mean)
ld2a$dist = "1000-10,000"


ld3 = ld[ld$dist > 10000 & ld$dist < 100000,]
ld3a = aggregate(ld3$R2,by=list(ld3$CHR_A), FUN=mean)
ld3a$dist = "10,000-100,000"


ld4 = ld[ld$dist > 100000 & ld$dist < 1000000,]
ld4a = aggregate(ld4$R2,by=list(ld4$CHR_A), FUN=mean)
ld4a$dist = "100,000-1,000,000"

lda = rbind(ld1a,ld2a,ld3a,ld4a)
lda$chr = paste0("chr",lda[,1])


df3 = merge(lda,df2,by = "chr")
df3$dist = factor(df3$dist, levels=c("<1000","1000-10,000","10,000-100,000","100,000-1,000,000"))



outpl = ggplot(df3, aes(x=pi, y=x)) + geom_point() + facet_wrap(~ dist) + theme_bw() + 
  theme(text = element_text(size = 12),  axis.text = element_text(colour = "black"), panel.grid.minor = element_blank() , panel.grid.major = element_blank()) +
    xlab("Mean nucleotide diversity on chromosome") + ylab(expression(paste("Mean LD (", "r"^2, ") for distance class")))


ggsave(
  "LD_pi.pdf",
  plot = outpl,
  device = "pdf",
  path = "/plas1/george.sandler/marchantia_popgen_newgenome/ld/",
  scale = 1,
  width = 35,
  height = 13,
  units = c("cm"),
  dpi = 500,
  limitsize = TRUE,
) 

     