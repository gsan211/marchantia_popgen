

library(data.table)
library(ggplot2)
library(dplyr)

df = NULL
for(i in c(1:5)){

  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex/hmmibd/slim_facsex_90",i,".hmm_fract.txt") 
  #file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/newsubpop/hmmibd/slim_newsubpop_95_05_",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  df = rbind(df,dfi)

}



#ggplot(df, aes(x=fract_sites_IBD, y = ..scaled..,fill = as.factor(run) )) + theme_bw() + geom_density(alpha=0.3)+ geom_vline(xintercept = 0.1735537)

ggplot(df, aes(x=fract_sites_IBD, y = ..scaled..,fill = as.factor(run) )) + theme_bw() + geom_density(alpha=0.3)

ggplot(df, aes(x=fract_sites_IBD, y = N_generation,colour = as.factor(run) )) + theme_bw() + geom_point(alpha=0.3)







######################

df = NULL
for(i in c(1:5)){

  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex/test/slim_facsex_999",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  dfi$asex = 0.999
  df = rbind(df,dfi)
  
  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex/test/slim_facsex_99",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  dfi$asex = 0.99
  df = rbind(df,dfi)


  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex/test/slim_facsex_90",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  dfi$asex = 0.90
  df = rbind(df,dfi)
  
  
  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex/test/slim_facsex_00",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  dfi$asex = 0
  df = rbind(df,dfi)
}

ggplot(df, aes(x=fract_sites_IBD, y = ..scaled..,fill = as.factor(run) )) + theme_bw() + geom_density(alpha=0.3) +  facet_wrap(~ asex)


ggplot(df, aes(x=fract_sites_IBD, y = N_generation,colour = as.factor(run) )) + theme_bw() + geom_point(alpha=0.3)


df2 = df %>% group_by(run, asex) %>% summarise_at(vars(fract_sites_IBD), funs(mean, sd)) 
df3 = df2 %>% group_by(asex) %>% summarise_at(vars(mean,sd), funs(mean))
  
pl1 = ggplot(df, aes(x=fract_sites_IBD, y = ..scaled..,fill = as.factor(run) )) + theme_bw() + geom_density(alpha=0.3) +  facet_wrap(~ asex) + xlab("Fraction of genome IBD per sample pair") + ylab("Scaled Density") +
theme(text = element_text(size = 18),axis.text=element_text(colour = "black"),legend.position = "none")

##################
ggsave(
  "slim_ibd_10k.pdf",
  plot = pl1,
  device = "pdf",
  path = "/plas1/george.sandler/marchantia_popgen_newgenome/slim/",
  scale = 1,
  width = 30,
  height = 25,
  units = c("cm"),
  dpi = 500,
  limitsize = TRUE,
) 
 
 
df2 = df %>% group_by(run, asex) %>% summarise_at(vars(fract_sites_IBD), funs(mean, sd)) 
df3 = df2 %>% group_by(asex) %>% summarise_at(vars(mean,sd), funs(mean))
  
  
##########################


df = NULL
for(i in c(1:5)){

  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex2k/test/slim_facsex2k_999",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  dfi$asex = 0.999
  df = rbind(df,dfi)
  
  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex2k/test/slim_facsex2k_99",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  dfi$asex = 0.99
  df = rbind(df,dfi)


  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex2k/test/slim_facsex2k_90",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  dfi$asex = 0.9
  df = rbind(df,dfi)
  
  
  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex2k/test/slim_facsex2k_00",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  dfi$asex = 0
  df = rbind(df,dfi)
}

pl2 = ggplot(df, aes(x=fract_sites_IBD, y = ..scaled..,fill = as.factor(run) )) + theme_bw() + geom_density(alpha=0.3) +  facet_wrap(~ asex) + xlab("Fraction of genome IBD per sample pair") + ylab("Scaled Density") +
theme(text = element_text(size = 18),axis.text=element_text(colour = "black"),legend.position = "none")



ggsave(
  "slim_ibd_2k.pdf",
  plot = pl2,
  device = "pdf",
  path = "/plas1/george.sandler/marchantia_popgen_newgenome/slim/",
  scale = 1,
  width = 30,
  height = 25,
  units = c("cm"),
  dpi = 500,
  limitsize = TRUE,
) 
 
 
df2 = df %>% group_by(run, asex) %>% summarise_at(vars(fract_sites_IBD), funs(mean, sd)) 
df3 = df2 %>% group_by(asex) %>% summarise_at(vars(mean,sd), funs(mean))
  


  
##############################

  #file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/admix_subpop/hmmibd/slim_subpop_ad50k_",i,".hmm_fract.txt")
  #file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/admix_subpop/hmmibd/slim_subpop_noad_",i,".hmm_fract.txt")
  #file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/admix_subpop/hmmibd/slim_subpop_noad_large_",i,".hmm_fract.txt")
  #file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/admix_subpop/hmmibd/slim_subpop_ad65k_10x_input",i,".hmm_fract.txt")