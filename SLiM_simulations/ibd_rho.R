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



df = NULL
for(i in c(1:5)){

  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/low_rho/hmmibd/slim_rec_09_",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  dfi$asex = 1e-09
  df = rbind(df,dfi)
  
  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/low_rho/hmmibd/slim_rec_10_",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  dfi$asex = 1e-10
  df = rbind(df,dfi)


  file_in = paste0("/plas1/george.sandler/marchantia_popgen_newgenome/slim/fac_sex/hmmibd/slim_facsex_00",i,".hmm_fract.txt") 
  dfi = fread(file_in, header = T)
  dfi$run = i
  dfi$asex = 1e-08
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
  "slim_ibd_rho.pdf",
  plot = pl1,
  device = "pdf",
  path = "/plas1/george.sandler/marchantia_popgen_newgenome/slim/",
  scale = 1,
  width = 43,
  height = 15,
  units = c("cm"),
  dpi = 500,
  limitsize = TRUE,
) 
 
 
df2 = df %>% group_by(run, asex) %>% summarise_at(vars(fract_sites_IBD), funs(mean, sd)) 
df3 = df2 %>% group_by(asex) %>% summarise_at(vars(mean,sd), funs(mean))
  